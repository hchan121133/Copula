############################################################
# cvine_pay.R  (METHOD = PAY)
# - Synthetic data generation included (raw_df, cum_df)
# - C-vine copula per (g,t) with hurdle margins (parametric + ECDF fallback)
# - Progress + parallel execution
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(fitdistrplus)
  library(VineCopula)
  library(future)
  library(future.apply)
  library(progressr)
})

options(scipen = 999)

############################################################
# [DATA GEN] 가상 데이터 생성 (a,b,c 동시지급 패턴 + 5억 초과 튜닝)
# a = y1~y6, b = y7, c = y8
# - 생성 결과: raw_df, cum_df (g 포함)
############################################################
if (!exists("raw_df") || !exists("cum_df")) {
  
  set.seed(20260206)
  
  N     <- 20000
  Tmax  <- 20
  LIMIT <- 5e8
  
  amt_y1 <- 1000 * 10^4
  amt_y2 <- 5000 * 10^4
  amt_y3 <- 10000 * 10^4
  amt_y4 <- 1000 * 10^4
  amt_y5 <- 2000 * 10^4
  amt_y6 <- 2000 * 10^4
  amt_y7 <- 10 * 10^4
  amt_y8 <- c(20*10^4, 40*10^4, 60*10^4, 100*10^4, 200*10^4, 600*10^4, 1000*10^4)
  
  CAP_N1_PER_YEAR <- 8
  logit <- function(x) 1/(1+exp(-x))
  
  id_tbl <- tibble::tibble(
    id        = sprintf("ID%06d", 1:N),
    sex       = sample(c("M","F"), N, TRUE),
    issue_age = sample(15:69, N, TRUE)
  ) %>%
    dplyr::mutate(
      age_band = cut(
        issue_age,
        breaks = c(0,19,29,39,49,59,69,120),
        labels = c("15-19","20-29","30-39","40-49","50-59","60-69","70+")
      )
    )
  
  latent_tbl <- tibble::tibble(
    id       = id_tbl$id,
    Z_common = rnorm(N, sd = 1.3),
    Z_surg   = rnorm(N, sd = 1.3)
  )
  
  panel <- id_tbl %>%
    tidyr::crossing(duration = 1:Tmax) %>%
    dplyr::mutate(attained_age = issue_age + duration - 1) %>%
    dplyr::left_join(latent_tbl, by = "id")
  
  tmp <- panel %>%
    dplyr::mutate(
      lambda1 = exp(-4.6 + 0.05*(attained_age-40) + 0.03*(duration-1) + 0.85*Z_common),
      n1 = pmin(rpois(dplyr::n(), lambda = lambda1), CAP_N1_PER_YEAR),
      y1 = n1 * amt_y1,
      
      n2 = rbinom(dplyr::n(), 1, prob = logit(-5.6 + 0.06*(attained_age-40) + 0.75*Z_common)),
      y2 = n2 * amt_y2,
      
      n3 = rbinom(dplyr::n(), 1, prob = logit(-6.0 + 0.07*(attained_age-40) + 0.95*Z_common)),
      y3 = n3 * amt_y3,
      
      n4 = rbinom(dplyr::n(), 1, prob = logit(-5.8 + 0.06*(attained_age-40) + 0.80*Z_common)),
      y4 = n4 * amt_y4,
      
      n5 = rbinom(dplyr::n(), 1, prob = logit(-5.6 + 0.07*(attained_age-45) + 0.25*(sex=="M") + 0.85*Z_common)),
      y5 = n5 * amt_y5,
      
      n6 = rbinom(dplyr::n(), 1, prob = logit(-5.9 + 0.07*(attained_age-45) + 0.30*(sex=="M") + 0.90*Z_common)),
      y6 = n6 * amt_y6,
      
      a_pay   = y1 + y2 + y3 + y4 + y5 + y6,
      a_event = as.integer(a_pay > 0)
    )
  
  tmp <- tmp %>%
    dplyr::mutate(
      p7_admit = logit(-3.4 + 0.04*(attained_age-40) + 0.65*Z_common + 1.60*a_event),
      admit7   = rbinom(dplyr::n(), 1, prob = p7_admit),
      
      lam_days = pmax(1, 6 + 6*a_event + 2.0*Z_common),
      n7 = ifelse(admit7==1, rpois(dplyr::n(), lambda = lam_days), 0),
      n7 = pmin(n7, 30),
      
      y7    = n7 * amt_y7,
      b_pay = y7
    )
  
  for (k in 1:7) {
    base_p <- 0.020 + 0.003*k
    high_class_boost <- ifelse(k >= 6, 0.60, 0.00)
    
    tmp[[paste0("p8_",k)]] <- logit(
      qlogis(base_p) +
        0.70*tmp$Z_surg +
        0.55*tmp$Z_common +
        1.40*tmp$a_event +
        high_class_boost
    )
    tmp[[paste0("n8_",k)]] <- rbinom(nrow(tmp), 1, prob = tmp[[paste0("p8_",k)]])
    tmp[[paste0("y8_",k)]] <- tmp[[paste0("n8_",k)]] * amt_y8[k]
  }
  
  tmp <- tmp %>%
    dplyr::mutate(
      y8      = y8_1 + y8_2 + y8_3 + y8_4 + y8_5 + y8_6 + y8_7,
      c_pay   = y8,
      y_total = a_pay + b_pay + c_pay
    )
  
  raw_df <- tmp %>%
    dplyr::arrange(id, duration) %>%
    dplyr::group_by(id) %>%
    dplyr::mutate(
      S_cum      = cumsum(y_total),
      exceed_5e8 = (S_cum > 5e8)
    ) %>%
    dplyr::ungroup() %>%
    # [SELECT-DATA-01]
    dplyr::select(
      id, sex, issue_age, age_band, duration, attained_age,
      n1,n2,n3,n4,n5,n6, y1,y2,y3,y4,y5,y6, a_pay,
      n7, y7, b_pay,
      dplyr::starts_with("n8_"), y8, c_pay,
      y_total, S_cum, exceed_5e8
    )
  
  exceed_id_rate <- raw_df %>%
    dplyr::group_by(id) %>%
    dplyr::summarise(exceed_any = any(S_cum > 5e8), .groups = "drop") %>%
    dplyr::summarise(rate = mean(exceed_any), n_exceed = sum(exceed_any))
  print(exceed_id_rate)
  
  check_tbl <- raw_df %>%
    dplyr::mutate(a_event = (a_pay > 0)) %>%
    dplyr::summarise(
      b_rate_when_a     = mean(b_pay > 0 & a_event),
      b_rate_when_not_a = mean(b_pay > 0 & !a_event),
      c_rate_when_a     = mean(c_pay > 0 & a_event),
      c_rate_when_not_a = mean(c_pay > 0 & !a_event)
    )
  print(check_tbl)
  
  cum_df <- raw_df %>%
    dplyr::arrange(id, duration) %>%
    dplyr::group_by(id) %>%
    dplyr::mutate(
      A_cum = cumsum(a_pay),
      B_cum = cumsum(b_pay),
      C_cum = cumsum(c_pay),
      S_cum_abc = A_cum + B_cum + C_cum
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(g = paste0(sex, "_", age_band))
}

############################################################
# Common functions (margins + C-vine + simulation + outputs)
############################################################

.make_ecdf_FQ <- function(x_pos) {
  x_pos <- sort(as.numeric(x_pos))
  x_pos <- x_pos[is.finite(x_pos) & x_pos > 0]
  if (length(x_pos) == 0) {
    return(list(
      F_pos = function(z) rep(0, length(z)),
      Q_pos = function(p) rep(0, length(p))
    ))
  }
  Fn <- stats::ecdf(x_pos)
  F_pos <- function(z) pmin(pmax(Fn(as.numeric(z)), 1e-10), 1 - 1e-10)
  Q_pos <- function(p) {
    p <- pmin(pmax(as.numeric(p), 1e-10), 1 - 1e-10)
    as.numeric(stats::quantile(x_pos, probs = p, type = 8, names = FALSE))
  }
  list(F_pos = F_pos, Q_pos = Q_pos)
}

.safe_fitdist <- function(x_pos, dist_name) {
  x_pos <- as.numeric(x_pos)
  x_pos <- x_pos[is.finite(x_pos) & x_pos > 1e-8]
  if (length(x_pos) < 5) return(NULL)
  if (!is.finite(sd(x_pos)) || sd(x_pos) < 1e-10) return(NULL)
  
  start_list <- NULL
  extra_args1 <- list()
  extra_args2 <- list()
  
  if (dist_name == "lnorm") {
    lx <- log(x_pos)
    if (!all(is.finite(lx))) return(NULL)
    sdl <- sd(lx); if (!is.finite(sdl) || sdl < 1e-6) sdl <- 1e-3
    start_list <- list(meanlog = mean(lx), sdlog = sdl)
    extra_args1$optim.method <- "Nelder-Mead"
    extra_args2$optim.method <- "BFGS"
    
  } else if (dist_name == "gamma") {
    m <- mean(x_pos); v <- var(x_pos)
    if (!is.finite(v) || v <= 0) return(NULL)
    shape0 <- max(0.2, m^2 / v)
    rate0  <- max(1e-6, m / v)
    start_list <- list(shape = shape0, rate = rate0)
    extra_args1$optim.method <- "L-BFGS-B"
    extra_args1$lower <- c(shape = 1e-8, rate = 1e-12)
    extra_args2$optim.method <- "Nelder-Mead"
    
  } else if (dist_name == "weibull") {
    m <- mean(x_pos)
    start_list <- list(shape = 1.2, scale = max(1e-6, m))
    extra_args1$optim.method <- "L-BFGS-B"
    extra_args1$lower <- c(shape = 1e-8, scale = 1e-12)
    extra_args2$optim.method <- "Nelder-Mead"
  } else {
    return(NULL)
  }
  
  .quiet_fit <- function(extra_args) {
    tryCatch(
      suppressWarnings(
        suppressMessages(
          do.call(
            fitdistrplus::fitdist,
            c(list(data = x_pos, distr = dist_name, start = start_list), extra_args)
          )
        )
      ),
      error = function(e) NULL
    )
  }
  
  fit1 <- .quiet_fit(extra_args1); if (!is.null(fit1)) return(fit1)
  fit2 <- .quiet_fit(extra_args2); if (!is.null(fit2)) return(fit2)
  NULL
}

fit_hurdle_margin_FQ <- function(x, candidates = c("lnorm","gamma","weibull")) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  p0 <- mean(x == 0)
  
  x_pos <- x[x > 0]
  x_pos <- x_pos[is.finite(x_pos) & x_pos > 1e-8]
  
  fits <- list()
  aic_tbl <- tibble(dist = character(), AIC = double())
  
  for (d in candidates) {
    fit_obj <- .safe_fitdist(x_pos, d)
    if (!is.null(fit_obj)) {
      fits[[d]] <- fit_obj
      aic_tbl <- dplyr::bind_rows(aic_tbl, tibble(dist = d, AIC = AIC(fit_obj)))
    }
  }
  
  if (nrow(aic_tbl) == 0) {
    best_dist <- NA_character_
    best_fit  <- NULL
  } else {
    best_dist <- aic_tbl %>% dplyr::arrange(AIC) %>% dplyr::slice(1) %>% dplyr::pull(dist)
    best_fit  <- fits[[best_dist]]
  }
  
  if (!is.null(best_fit)) {
    pars <- best_fit$estimate
    F_pos <- switch(
      best_dist,
      "lnorm"   = function(z) plnorm(z, meanlog = pars["meanlog"], sdlog = pars["sdlog"]),
      "gamma"   = function(z) {
        rate_val <- if ("rate" %in% names(pars)) pars["rate"] else 1/pars["scale"]
        pgamma(z, shape = pars["shape"], rate = rate_val)
      },
      "weibull" = function(z) pweibull(z, shape = pars["shape"], scale = pars["scale"]),
      .make_ecdf_FQ(x_pos)$F_pos
    )
    Q_pos <- switch(
      best_dist,
      "lnorm"   = function(p) qlnorm(p, meanlog = pars["meanlog"], sdlog = pars["sdlog"]),
      "gamma"   = function(p) {
        rate_val <- if ("rate" %in% names(pars)) pars["rate"] else 1/pars["scale"]
        qgamma(p, shape = pars["shape"], rate = rate_val)
      },
      "weibull" = function(p) qweibull(p, shape = pars["shape"], scale = pars["scale"]),
      .make_ecdf_FQ(x_pos)$Q_pos
    )
    note <- "parametric"
  } else {
    ec <- .make_ecdf_FQ(x_pos)
    F_pos <- ec$F_pos
    Q_pos <- ec$Q_pos
    note <- "ecdf_fallback"
  }
  
  F_fun <- function(z) {
    z <- as.numeric(z)
    out <- rep(NA_real_, length(z))
    out[z < 0]  <- 0
    out[z == 0] <- p0
    idx <- which(z > 0)
    out[idx] <- p0 + (1 - p0) * F_pos(z[idx])
    pmin(pmax(out, 1e-10), 1 - 1e-10)
  }
  
  Q_fun <- function(u) {
    u <- pmin(pmax(as.numeric(u), 1e-10), 1 - 1e-10)
    out <- rep(0, length(u))
    idx <- which(u > p0)
    if (length(idx) > 0) {
      p <- (u[idx] - p0) / max(1e-12, (1 - p0))
      out[idx] <- Q_pos(p)
    }
    out
  }
  
  list(p0=p0, best_dist=best_dist, best_fit=best_fit, aic_tbl=aic_tbl, note=note, F_fun=F_fun, Q_fun=Q_fun)
}

get_nearest_t_data <- function(df, g_value, t_value, t_col = "duration") {
  d0 <- df %>% dplyr::filter(g == g_value, .data[[t_col]] == t_value)
  if (nrow(d0) > 0) return(d0)
  ts <- df %>% dplyr::filter(g == g_value) %>% dplyr::distinct(.data[[t_col]]) %>% dplyr::pull()
  if (length(ts) == 0) return(d0)
  t_near <- ts[which.min(abs(ts - t_value))]
  df %>% dplyr::filter(g == g_value, .data[[t_col]] == t_near)
}

fit_cvine_3d <- function(U, root = 1, familyset = NA) {
  U <- as.matrix(U)
  U <- pmin(pmax(U, 1e-10), 1 - 1e-10)
  
  ord <- c(root, setdiff(1:3, root))
  Uo <- U[, ord, drop = FALSE]
  
  rvm <- tryCatch(
    suppressWarnings(
      suppressMessages(
        VineCopula::RVineStructureSelect(
          Uo,
          familyset = familyset,
          type = 1,
          selectioncrit = "AIC",
          indeptest = TRUE,
          level = 0.05
        )
      )
    ),
    error = function(e) NULL
  )
  
  if (is.null(rvm)) return(list(rvm=NULL, ord=ord, indep=TRUE))
  list(rvm=rvm, ord=ord, indep=FALSE)
}

sim_cvine_3d <- function(cvine_fit, n) {
  if (!is.null(cvine_fit$indep) && cvine_fit$indep) {
    U <- matrix(runif(n * 3), n, 3)
    return(pmin(pmax(U, 1e-10), 1 - 1e-10))
  }
  Uo_sim <- VineCopula::RVineSim(n, cvine_fit$rvm)
  inv_ord <- match(1:3, cvine_fit$ord)
  U_sim <- Uo_sim[, inv_ord, drop = FALSE]
  pmin(pmax(U_sim, 1e-10), 1 - 1e-10)
}

wallet_exhaust_year <- function(S_cum_vec, limit) {
  idx <- which(S_cum_vec > limit)
  if (length(idx) == 0) return(NA_integer_)
  min(idx)
}

risk_summary <- function(S_mat, limit) {
  Tmax <- ncol(S_mat)
  S_T <- S_mat[, Tmax]
  ex_year <- apply(S_mat, 1, wallet_exhaust_year, limit = limit)
  q95 <- as.numeric(quantile(S_T, 0.95, type = 8))
  tibble(
    prob_ex_by_T = mean(!is.na(ex_year)),
    mean_S_T = mean(S_T),
    VaR95_S_T = q95,
    TVaR95_S_T = mean(S_T[S_T >= q95])
  )
}

exhaust_curve <- function(S_mat, limit) {
  Tmax <- ncol(S_mat)
  ex_year <- apply(S_mat, 1, wallet_exhaust_year, limit = limit)
  prob_by <- sapply(1:Tmax, function(t) mean(!is.na(ex_year) & ex_year <= t))
  prob_at <- sapply(1:Tmax, function(t) mean(ex_year == t))
  tibble(duration = 1:Tmax, prob_exhaust_by_t = prob_by, prob_exhaust_at_t = prob_at)
}

model_to_row_cvine <- function(method, g, t, model) {
  cop_note <- if (!is.null(model$cop$indep) && model$cop$indep) "indep_fallback" else "cvine"
  tibble(
    method = method, g = g, duration = t, nobs = model$nobs,
    copula_family = "C-vine", cop_note = cop_note,
    cop_AIC = if (cop_note == "cvine") model$cop$rvm$aic else NA_real_,
    cop_logLik = if (cop_note == "cvine") model$cop$rvm$logLik else NA_real_,
    A_note = model$mA$note, A_p0 = model$mA$p0, A_best = model$mA$best_dist,
    B_note = model$mB$note, B_p0 = model$mB$p0, B_best = model$mB$best_dist,
    C_note = model$mC$note, C_p0 = model$mC$p0, C_best = model$mC$best_dist
  )
}

extract_cvine_pairs <- function(rvm, ord, method, g, t) {
  if (is.null(rvm)) {
    return(tibble(
      method = method, g = g, duration = t,
      pair = NA_character_, cond = NA_character_,
      family = NA_integer_, family_name = NA_character_,
      par = NA_real_, par2 = NA_real_,
      note = "indep_fallback"
    ))
  }
  fam_mat  <- rvm$family
  par_mat  <- rvm$par
  par2_mat <- rvm$par2
  fam_name <- function(f) ifelse(is.na(f), NA_character_, VineCopula::BiCopName(f))
  
  root <- ord[1]; o2 <- ord[2]; o3 <- ord[3]
  
  tibble(
    method = method, g = g, duration = t,
    pair = c(paste0(root,"-",o2), paste0(root,"-",o3), paste0(o2,"-",o3)),
    cond = c("none","none",paste0("|",root)),
    family = c(fam_mat[2,1], fam_mat[3,1], fam_mat[3,2]),
    par = c(par_mat[2,1], par_mat[3,1], par_mat[3,2]),
    par2 = c(par2_mat[2,1], par2_mat[3,1], par2_mat[3,2]),
    family_name = c(fam_name(fam_mat[2,1]), fam_name(fam_mat[3,1]), fam_name(fam_mat[3,2])),
    note = "cvine"
  ) %>%
    dplyr::mutate(
      pair = dplyr::case_when(
        pair %in% c("1-2","2-1") ~ "A-B",
        pair %in% c("1-3","3-1") ~ "A-C",
        pair %in% c("2-3","3-2") ~ "B-C",
        TRUE ~ pair
      ),
      cond = dplyr::case_when(
        cond=="|1" ~ "|A",
        cond=="|2" ~ "|B",
        cond=="|3" ~ "|C",
        TRUE ~ cond
      )
    )
}

############################################################
# METHOD=PAY 데이터 구성  (연간지급액)
############################################################
method <- "pay"
method_dat <- raw_df %>%
  dplyr::mutate(g = paste0(sex, "_", age_band)) %>%
  # [SELECT-01]
  dplyr::select(id, g, duration, A = a_pay, B = b_pay, C = c_pay)

.fit_models_generic <- function(method_dat, g_value, Tmax, candidates, cvine_root = 1, familyset = NA) {
  models <- vector("list", Tmax)
  dat0_all <- method_dat %>% dplyr::filter(g == g_value)
  
  for (t in 1:Tmax) {
    dat <- get_nearest_t_data(dat0_all, g_value, t, "duration")
    
    if (nrow(dat) == 0) {
      mA <- fit_hurdle_margin_FQ(rep(0, 10), candidates)
      mB <- fit_hurdle_margin_FQ(rep(0, 10), candidates)
      mC <- fit_hurdle_margin_FQ(rep(0, 10), candidates)
      cop <- list(rvm = NULL, ord = c(1,2,3), indep = TRUE)
      models[[t]] <- list(mA=mA, mB=mB, mC=mC, cop=cop, nobs=0)
      next
    }
    
    mA <- fit_hurdle_margin_FQ(dat$A, candidates)
    mB <- fit_hurdle_margin_FQ(dat$B, candidates)
    mC <- fit_hurdle_margin_FQ(dat$C, candidates)
    
    U <- cbind(mA$F_fun(dat$A), mB$F_fun(dat$B), mC$F_fun(dat$C))
    cop <- fit_cvine_3d(U, root = cvine_root, familyset = familyset)
    
    models[[t]] <- list(mA=mA, mB=mB, mC=mC, cop=cop, nobs=nrow(dat))
  }
  models
}

.simulate_generic <- function(models, n_paths, Tmax) {
  
  A_mat <- matrix(0, n_paths, Tmax)
  B_mat <- matrix(0, n_paths, Tmax)
  C_mat <- matrix(0, n_paths, Tmax)
  
  for (t in 1:Tmax) {
    U_sim <- sim_cvine_3d(models[[t]]$cop, n_paths)
    A_mat[, t] <- pmax(0, models[[t]]$mA$Q_fun(U_sim[, 1]))
    B_mat[, t] <- pmax(0, models[[t]]$mB$Q_fun(U_sim[, 2]))
    C_mat[, t] <- pmax(0, models[[t]]$mC$Q_fun(U_sim[, 3]))
  }
  
  A_cum <- t(apply(A_mat, 1, cumsum))
  B_cum <- t(apply(B_mat, 1, cumsum))
  C_cum <- t(apply(C_mat, 1, cumsum))
  S_cum <- A_cum + B_cum + C_cum
  
  list(A=A_cum, B=B_cum, C=C_cum, S=S_cum)
}

run_method <- function(method_dat,
                       method_name = "pay",
                       candidates = c("lnorm","gamma","weibull"),
                       Tmax = max(method_dat$duration),
                       n_paths = 5000,      # ★ 시뮬레이션 횟수(여기)
                       limit = 5e8,
                       cvine_root = 1,
                       familyset = NA,
                       workers = NULL) {
  
  g_list <- sort(unique(method_dat$g))
  n_tasks <- length(g_list)
  
  if (is.null(workers)) {
    workers <- max(1, min(n_tasks, parallel::detectCores(logical = FALSE) - 1))
  } else {
    workers <- max(1, min(workers, n_tasks))
  }
  
  message("=== ", method_name, " / C-vine ===")
  message("groups=", n_tasks, ", workers=", workers, ", Tmax=", Tmax, ", n_paths=", n_paths)
  
  future::plan(future::multisession, workers = workers)
  progressr::handlers(global = TRUE)
  progressr::handlers("txtprogressbar")
  
  one_g <- function(g) {
    t0 <- Sys.time()
    models <- .fit_models_generic(method_dat, g, Tmax, candidates, cvine_root, familyset)
    sim <- .simulate_generic(models, n_paths, Tmax)
    
    model_tbl <- dplyr::bind_rows(lapply(1:Tmax, function(t) model_to_row_cvine(method_name, g, t, models[[t]])))
    pair_tbl  <- dplyr::bind_rows(lapply(1:Tmax, function(t) {
      cop <- models[[t]]$cop
      rvm <- if (!is.null(cop$indep) && cop$indep) NULL else cop$rvm
      extract_cvine_pairs(rvm, cop$ord, method_name, g, t)
    }))
    curve_tbl   <- exhaust_curve(sim$S, limit) %>% dplyr::mutate(method = method_name, g = g)
    summary_tbl <- risk_summary(sim$S, limit) %>% dplyr::mutate(method = method_name, g = g)
    
    list(
      g=g,
      model_tbl=model_tbl,
      pair_tbl=pair_tbl,
      curve_tbl=curve_tbl,
      summary_tbl=summary_tbl,
      elapsed_sec = as.numeric(difftime(Sys.time(), t0, units = "secs"))
    )
  }
  
  t0_all <- Sys.time()
  out_list <- progressr::with_progress({
    p <- progressr::progressor(steps = n_tasks)
    future.apply::future_lapply(
      g_list,
      function(g) { res <- one_g(g); p(g); res },
      future.seed = TRUE
    )
  })
  
  model_tbl   <- dplyr::bind_rows(lapply(out_list, `[[`, "model_tbl"))
  pair_tbl    <- dplyr::bind_rows(lapply(out_list, `[[`, "pair_tbl"))
  curve_tbl   <- dplyr::bind_rows(lapply(out_list, `[[`, "curve_tbl"))
  summary_tbl <- dplyr::bind_rows(lapply(out_list, `[[`, "summary_tbl"))
  
  fallback_tbl <- model_tbl %>%
    dplyr::group_by(method, g) %>%
    dplyr::summarise(
      A_ecdf_rate = mean(A_note == "ecdf_fallback"),
      B_ecdf_rate = mean(B_note == "ecdf_fallback"),
      C_ecdf_rate = mean(C_note == "ecdf_fallback"),
      indep_cvine_rate = mean(cop_note == "indep_fallback"),
      .groups = "drop"
    )
  
  time_tbl <- tibble(
    method = method_name,
    g = vapply(out_list, `[[`, "", "g"),
    elapsed_sec = vapply(out_list, `[[`, 0.0, "elapsed_sec")
  ) %>% dplyr::arrange(g)
  
  message("DONE ", method_name, ". total(sec)=", round(as.numeric(difftime(Sys.time(), t0_all, units = "secs")), 1))
  
  list(
    model_tbl=model_tbl, pair_tbl=pair_tbl, curve_tbl=curve_tbl, summary_tbl=summary_tbl,
    fallback_tbl=fallback_tbl, time_tbl=time_tbl,
    meta=list(method=method_name, Tmax=Tmax, n_paths=n_paths, limit=limit, cvine_root=cvine_root, familyset=familyset, workers=workers)
  )
}

############################################################
# 실행 + 저장
############################################################
Tmax_run <- max(method_dat$duration)
LIMIT_run <- 5e8
N_PATHS <- 5000        # ★ 시뮬레이션 횟수(여기)
ROOT <- 1

out_pay <- run_method(method_dat, method_name="pay", Tmax=Tmax_run, n_paths=N_PATHS, limit=LIMIT_run, cvine_root=ROOT)

saveRDS(out_pay, "out_pay.rds")
write.csv(out_pay$model_tbl,   "pay_model_tbl.csv",   row.names = FALSE)
write.csv(out_pay$pair_tbl,    "pay_pair_tbl.csv",    row.names = FALSE)
write.csv(out_pay$curve_tbl,   "pay_curve_tbl.csv",   row.names = FALSE)
write.csv(out_pay$summary_tbl, "pay_summary_tbl.csv", row.names = FALSE)
write.csv(out_pay$fallback_tbl,"pay_fallback_tbl.csv",row.names = FALSE)
write.csv(out_pay$time_tbl,    "pay_time_tbl.csv",    row.names = FALSE)

message("Saved: out_pay.rds + pay_*.csv")
