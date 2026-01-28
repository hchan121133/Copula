############################################################

# FINAL B) Time-path (3*T=60D) "C-VINE" + Simulation

#

# - For each group g:

#   1) Build 60D vector per id: A_inc_1..T, B_inc_1..T, C_inc_1..T

#   2) Fit hurdle margins for each dimension

#   3) Fit C-vine copula on 60D U (forced)

#   4) Simulate 60D paths -> cumulate -> exceed probability

# - Export:

#   (1) exceed summary CSV

#   (2) exceed-by-year CSV

#   (3) margin fit summary CSV

#   (4) vine family counts CSV

# - Includes "How to confirm it's C-vine" section

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



# =========================================================

# [B-0] USER CONTROLS

# =========================================================

LIMIT <- 10e8

Tmax  <- 20



min_pos_n <- 80



n_sim <- 1000         # ★ 시뮬레이션 횟수: 여기서 지정 ★



TRUNC_LEVEL <- 4        # ★ 60차원 계산량 제어 레버 ★ (2~3 권장)

familyset   <- c(1, 3, 4, 5, 6, 7,8,13, 14, 16)

out_dir2 <- "result_Model2"



# =========================================================

# [B-1] INPUT CHECK

# =========================================================

stopifnot(exists("cum_df"))

need_cols <- c("id","g","duration","A_cum","B_cum","C_cum")

stopifnot(all(need_cols %in% names(cum_df)))



if (!dir.exists(out_dir2)) dir.create(out_dir2, recursive = TRUE)

stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")



# =========================================================

# [B-2] CUM -> INC (long)

# =========================================================

inc_long <- cum_df %>%
  
  dplyr::arrange(id, duration) %>%
  
  dplyr::group_by(id) %>%
  
  dplyr::mutate(
    
    A_inc = A_cum - dplyr::lag(A_cum, default = 0),
    
    B_inc = B_cum - dplyr::lag(B_cum, default = 0),
    
    C_inc = C_cum - dplyr::lag(C_cum, default = 0)
    
  ) %>%
  
  dplyr::ungroup() %>%
  
  dplyr::mutate(
    
    A_inc = pmax(A_inc, 0),
    
    B_inc = pmax(B_inc, 0),
    
    C_inc = pmax(C_inc, 0)
    
  ) %>%
  
  dplyr::select(id, g, duration, A_inc, B_inc, C_inc)



# =========================================================

# [B-3] HURDLE UTILS (A와 동일)

# =========================================================

fit_hurdle_margin <- function(x, min_pos = 80) {
  
  x <- x[is.finite(x)]
  
  if (length(x) == 0) stop("empty x")
  
  
  
  p0 <- mean(x <= 0)
  
  y  <- x[x > 0]
  
  
  
  if (length(y) < min_pos || sd(y) == 0) {
    
    return(list(p0 = p0, pos_model = list(type = "ecdf", y = sort(y))))
    
  }
  
  
  
  cand <- list()
  
  cand[["lnorm"]] <- tryCatch({
    
    fit <- fitdist(y, "lnorm")
    
    list(type="lnorm", fit=fit, aic=AIC(fit))
    
  }, error=function(e) NULL)
  
  
  
  cand[["gamma"]] <- tryCatch({
    
    fit <- fitdist(y, "gamma")
    
    list(type="gamma", fit=fit, aic=AIC(fit))
    
  }, error=function(e) NULL)
  
  
  
  cand[["weibull"]] <- tryCatch({
    
    fit <- fitdist(y, "weibull")
    
    list(type="weibull", fit=fit, aic=AIC(fit))
    
  }, error=function(e) NULL)
  
  
  
  cand <- cand[!vapply(cand, is.null, logical(1))]
  
  if (length(cand) == 0) {
    
    return(list(p0 = p0, pos_model = list(type = "ecdf", y = sort(y))))
    
  }
  
  
  
  best <- cand[[ which.min(vapply(cand, `[[`, numeric(1), "aic")) ]]
  
  list(p0 = p0, pos_model = best)
  
}



cdf_pos <- function(z, pos_model) {
  
  if (pos_model$type == "ecdf") {
    
    y <- pos_model$y
    
    if (length(y) == 0) return(rep(0, length(z)))
    
    return(approx(x=y, y=seq_along(y)/length(y), xout=z, rule=2, ties="ordered")$y)
    
  }
  
  est <- pos_model$fit$estimate
  
  if (pos_model$type == "lnorm")   return(plnorm(z, meanlog = est["meanlog"], sdlog = est["sdlog"]))
  
  if (pos_model$type == "gamma")   return(pgamma(z, shape = est["shape"], rate = est["rate"]))
  
  if (pos_model$type == "weibull") return(pweibull(z, shape = est["shape"], scale = est["scale"]))
  
  stop("unknown pos_model type")
  
}



q_pos <- function(u, pos_model) {
  
  u <- pmin(pmax(u, 0), 1)
  
  if (pos_model$type == "ecdf") {
    
    y <- pos_model$y
    
    if (length(y) == 0) return(rep(0, length(u)))
    
    return(as.numeric(stats::quantile(y, probs = u, type = 8, names = FALSE)))
    
  }
  
  est <- pos_model$fit$estimate
  
  if (pos_model$type == "lnorm")   return(qlnorm(u, meanlog = est["meanlog"], sdlog = est["sdlog"]))
  
  if (pos_model$type == "gamma")   return(qgamma(u, shape = est["shape"], rate = est["rate"]))
  
  if (pos_model$type == "weibull") return(qweibull(u, shape = est["shape"], scale = est["scale"]))
  
  stop("unknown pos_model type")
  
}



pit_hurdle <- function(x, m) {
  
  p0 <- m$p0
  
  pos <- m$pos_model
  
  u <- numeric(length(x))
  
  is0 <- (x <= 0) | !is.finite(x)
  
  
  
  if (any(is0)) u[is0] <- runif(sum(is0), 0, max(p0, 1e-12))
  
  if (any(!is0)) {
    
    Fx <- cdf_pos(x[!is0], pos)
    
    Fx <- pmin(pmax(Fx, 0), 1)
    
    u[!is0] <- p0 + (1 - p0) * Fx
    
    u[!is0] <- pmin(pmax(u[!is0], 1e-10), 1 - 1e-10)
    
  }
  
  u
  
}



q_hurdle <- function(u, m) {
  
  p0 <- m$p0
  
  pos <- m$pos_model
  
  u <- pmin(pmax(u, 0), 1)
  
  
  
  x <- numeric(length(u))
  
  is0 <- (u <= p0)
  
  x[is0] <- 0
  
  
  
  if (any(!is0)) {
    
    uu <- (u[!is0] - p0) / (1 - p0)
    
    uu <- pmin(pmax(uu, 1e-10), 1 - 1e-10)
    
    x[!is0] <- q_pos(uu, pos)
    
  }
  
  x
  
}



# =========================================================

# [B-4] BUILD 60D WIDE TABLE PER GROUP

# =========================================================

make_wide_inc <- function(df_g_long, Tmax) {
  
  
  
  wideA <- df_g_long %>%
    
    dplyr::select(id, duration, A_inc) %>%
    
    tidyr::pivot_wider(names_from = duration, values_from = A_inc, values_fill = 0)
  
  names(wideA)[-1] <- paste0("A_inc_", names(wideA)[-1])
  
  
  
  wideB <- df_g_long %>%
    
    dplyr::select(id, duration, B_inc) %>%
    
    tidyr::pivot_wider(names_from = duration, values_from = B_inc, values_fill = 0)
  
  names(wideB)[-1] <- paste0("B_inc_", names(wideB)[-1])
  
  
  
  wideC <- df_g_long %>%
    
    dplyr::select(id, duration, C_inc) %>%
    
    tidyr::pivot_wider(names_from = duration, values_from = C_inc, values_fill = 0)
  
  names(wideC)[-1] <- paste0("C_inc_", names(wideC)[-1])
  
  
  
  out <- wideA %>%
    
    dplyr::left_join(wideB, by = "id") %>%
    
    dplyr::left_join(wideC, by = "id")
  
  
  
  wantA <- paste0("A_inc_", 1:Tmax)
  
  wantB <- paste0("B_inc_", 1:Tmax)
  
  wantC <- paste0("C_inc_", 1:Tmax)
  
  
  
  for (nm in c(wantA, wantB, wantC)) {
    
    if (!nm %in% names(out)) out[[nm]] <- 0
    
  }
  
  
  
  out %>% dplyr::select(id, dplyr::all_of(wantA), dplyr::all_of(wantB), dplyr::all_of(wantC))
  
}



# =========================================================

# [B-5] FIT 60D "FORCED C-VINE" PER GROUP

# - Primary: CVineStructureSelect(U, ...)

# - Fallback: RVineStructureSelect(U, type="CVine", trunclevel=TRUNC_LEVEL, ...)

#   (버전 이슈 대비: 실패 시 stop)

# =========================================================

fit_timevine_group_cvine <- function(df_wide, min_pos_n, familyset, trunc_level) {
  
  
  
  X <- df_wide %>% dplyr::select(-id)
  
  X <- as.data.frame(X)
  
  
  
  # (1) 각 차원별 hurdle 마진 적합
  
  margins <- vector("list", ncol(X))
  
  names(margins) <- names(X)
  
  for (j in seq_len(ncol(X))) {
    
    margins[[j]] <- fit_hurdle_margin(X[[j]], min_pos = min_pos_n)
    
  }
  
  
  
  # (2) PIT -> U
  
  U <- matrix(NA_real_, nrow=nrow(X), ncol=ncol(X))
  
  colnames(U) <- names(X)
  
  for (j in seq_len(ncol(X))) {
    
    U[,j] <- pit_hurdle(X[[j]], margins[[j]])
    
  }
  
  
  
  # (3) C-vine 강제 적합 시도
  
  vine <- NULL
  
  
  
  vine <- tryCatch({
    
    # 일부 버전에서는 trunclevel 인자가 없을 수 있음 → tryCatch로 보호
    
    CVineStructureSelect(
      
      U,
      
      familyset = familyset,
      
      selectioncrit = "AIC",
      
      indeptest = TRUE,
      
      method = "mle"
      
    )
    
  }, error = function(e) NULL)
  
  
  
  # (4) CVineStructureSelect가 실패하면, RVineStructureSelect + type="CVine"로 재시도
  
  if (is.null(vine)) {
    
    vine <- tryCatch({
      
      RVineStructureSelect(
        
        U,
        
        familyset = familyset,
        
        selectioncrit = "AIC",
        
        indeptest = TRUE,
        
        method = "mle",
        
        trunclevel = trunc_level,
        
        type = "CVine"
        
      )
      
    }, error = function(e) {
      
      stop("C-vine fit failed in both CVineStructureSelect and RVineStructureSelect(type='CVine').\n",
           
           "Please check VineCopula version / arguments.\n",
           
           "Original error: ", e$message)
      
    })
    
  }
  
  
  
  list(margins = margins, vine = vine, colnames = colnames(U))
  
}



# =========================================================

# [B-6] SIMULATE 60D PATHS

# =========================================================

simulate_timevine_group_B <- function(gv, mdl, n_sim, Tmax, LIMIT, seed) {
  
  set.seed(seed)
  
  
  
  U_sim <- RVineSim(n_sim, mdl$vine)
  
  colnames(U_sim) <- mdl$colnames
  
  
  
  X_sim <- matrix(0, nrow=n_sim, ncol=ncol(U_sim))
  
  colnames(X_sim) <- mdl$colnames
  
  
  
  for (j in seq_len(ncol(U_sim))) {
    
    X_sim[,j] <- pmax(q_hurdle(U_sim[,j], mdl$margins[[j]]), 0)
    
  }
  
  
  
  A_cols <- paste0("A_inc_", 1:Tmax)
  
  B_cols <- paste0("B_inc_", 1:Tmax)
  
  C_cols <- paste0("C_inc_", 1:Tmax)
  
  
  
  A_inc <- X_sim[, A_cols, drop=FALSE]
  
  B_inc <- X_sim[, B_cols, drop=FALSE]
  
  C_inc <- X_sim[, C_cols, drop=FALSE]
  
  
  
  A_cum <- t(apply(A_inc, 1, cumsum))
  
  B_cum <- t(apply(B_inc, 1, cumsum))
  
  C_cum <- t(apply(C_inc, 1, cumsum))
  
  S_cum <- A_cum + B_cum + C_cum
  
  
  
  exceed_any  <- apply(S_cum > LIMIT, 1, any)
  
  exceed_by_t <- colMeans(S_cum > LIMIT)
  
  
  
  list(
    
    g = gv,
    
    n_sim = n_sim,
    
    exceed_prob_any = mean(exceed_any),
    
    exceed_prob_by_t = exceed_by_t,
    
    S_cum_last = S_cum[,Tmax]
    
  )
  
}



# =========================================================

# [B-7] RUN FIT + SIM PER GROUP

# =========================================================

groups <- sort(unique(inc_long$g))



handlers(global = TRUE)

plan(multisession, workers = 2)



model_map_B <- NULL

sim_res_B <- NULL



with_progress({
  
  p <- progressor(along = groups)
  
  
  
  model_map_B <- future_lapply(groups, function(gv) {
    
    df_g <- inc_long %>% dplyr::filter(g == gv)
    
    df_wide <- make_wide_inc(df_g, Tmax = Tmax)
    
    
    
    mdl <- fit_timevine_group_cvine(
      
      df_wide = df_wide,
      
      min_pos_n = min_pos_n,
      
      familyset = familyset,
      
      trunc_level = TRUNC_LEVEL
      
    )
    
    
    
    p(sprintf("fit B g=%s", gv))
    
    mdl
    
  })
  
  
  
  names(model_map_B) <- groups
  
})



with_progress({
  
  p <- progressor(along = groups)
  
  
  
  sim_res_B <- future_lapply(groups, function(gv) {
    
    out <- simulate_timevine_group_B(
      
      gv = gv,
      
      mdl = model_map_B[[gv]],
      
      n_sim = n_sim,
      
      Tmax = Tmax,
      
      LIMIT = LIMIT,
      
      seed = 5000 + match(gv, groups)
      
    )
    
    p(sprintf("sim B g=%s", gv))
    
    out
    
  })
  
})



# =========================================================

# [B-8] OUTPUT: EXCEED RESULTS

# =========================================================

B_exceed_summary <- dplyr::bind_rows(lapply(sim_res_B, function(x) {
  
  tibble(
    
    method = "B",
    
    g = x$g,
    
    n_sim = x$n_sim,
    
    exceed_prob_any = x$exceed_prob_any,
    
    S_cum_last_p50 = as.numeric(quantile(x$S_cum_last, 0.50)),
    
    S_cum_last_p90 = as.numeric(quantile(x$S_cum_last, 0.90)),
    
    S_cum_last_p99 = as.numeric(quantile(x$S_cum_last, 0.99))
    
  )
  
})) %>% dplyr::arrange(dplyr::desc(exceed_prob_any))



B_exceed_by_t <- dplyr::bind_rows(lapply(sim_res_B, function(x) {
  
  tibble(
    
    method = "B",
    
    g = x$g,
    
    duration = 1:Tmax,
    
    exceed_prob_t = x$exceed_prob_by_t
    
  )
  
}))



B_exceed_summary_path <- file.path(out_dir2, paste0("B_exceed_summary_", stamp, ".csv"))

B_exceed_by_t_path    <- file.path(out_dir2, paste0("B_exceed_by_t_", stamp, ".csv"))



write.csv(B_exceed_summary, B_exceed_summary_path, row.names = FALSE)

write.csv(B_exceed_by_t,    B_exceed_by_t_path,    row.names = FALSE)



cat("Saved(B) exceed:\n", B_exceed_summary_path, "\n", B_exceed_by_t_path, "\n")



# =========================================================

# [B-9] OUTPUT: "MARGIN FIT SUMMARY" CSV (60D)

# =========================================================

extract_B_margin_summary <- function(model_map_B) {
  
  all_g <- names(model_map_B)
  
  out_list <- list()
  
  idx <- 1
  
  
  
  for (g in all_g) {
    
    mdl <- model_map_B[[g]]
    
    mlist <- mdl$margins
    
    vars <- names(mlist)
    
    
    
    for (v in vars) {
      
      m <- mlist[[v]]
      
      pos_type <- m$pos_model$type
      
      
      
      params <- NA_character_
      
      if (!is.null(m$pos_model$fit)) {
        
        est <- m$pos_model$fit$estimate
        
        params <- paste(names(est), round(as.numeric(est), 6), collapse="; ")
        
      }
      
      
      
      out_list[[idx]] <- data.frame(
        
        method = "B",
        
        g = g,
        
        var = v,
        
        p0 = m$p0,
        
        pos_type = pos_type,
        
        params = params,
        
        stringsAsFactors = FALSE
        
      )
      
      idx <- idx + 1
      
    }
    
  }
  
  
  
  dplyr::bind_rows(out_list) %>% dplyr::arrange(g, var)
  
}



B_margin_tbl <- extract_B_margin_summary(model_map_B)

B_margin_path <- file.path(out_dir2, paste0("B_margin_fit_summary_", stamp, ".csv"))

write.csv(B_margin_tbl, B_margin_path, row.names = FALSE)

cat("Saved(B) margin fit summary:\n", B_margin_path, "\n")



# =========================================================

# [B-10] OUTPUT: "VINE FAMILY COUNTS" CSV

# =========================================================

family_name <- function(code){
  
  nm <- VineCopula::BiCopName(code)
  
  ifelse(is.na(nm), as.character(code), nm)
  
}



extract_B_family_counts <- function(model_map_B) {
  
  all_g <- names(model_map_B)
  
  out_list <- list()
  
  idx <- 1
  
  
  
  for (g in all_g) {
    
    vine <- model_map_B[[g]]$vine
    
    
    
    fam_vec <- as.vector(vine$family)
    
    fam_vec <- fam_vec[fam_vec != 0]
    
    if (length(fam_vec) == 0) next
    
    
    
    tab <- sort(table(fam_vec), decreasing = TRUE)
    
    
    
    out_list[[idx]] <- data.frame(
      
      method = "B",
      
      g = g,
      
      family_code = as.integer(names(tab)),
      
      family_name = family_name(as.integer(names(tab))),
      
      n_edges = as.integer(tab),
      
      trunclevel = TRUNC_LEVEL,
      
      stringsAsFactors = FALSE
      
    )
    
    idx <- idx + 1
    
  }
  
  
  
  dplyr::bind_rows(out_list) %>%
    
    dplyr::arrange(g, dplyr::desc(n_edges))
  
}



B_family_tbl <- extract_B_family_counts(model_map_B)

B_family_path <- file.path(out_dir2, paste0("B_vine_family_counts_", stamp, ".csv"))

write.csv(B_family_tbl, B_family_path, row.names = FALSE)

cat("Saved(B) vine family counts:\n", B_family_path, "\n")



# =========================================================

# [B-11] "IS IT REALLY C-VINE?" 확인 방법

# =========================================================

# - 우리가 C-vine 강제 적합을 시도했고, 실패하면 중단(stop)하도록 해놨습니다.

# - 그래도 구조 확인하려면:

#   a) vine$Matrix 출력

#   b) summary(vine)

#   c) plot(vine)에서 Tree1이 star 형태(루트 1개가 다 연결)이면 C-vine 직관 확인



check_B_cvine <- function(model_map_B, g) {
  
  vine <- model_map_B[[g]]$vine
  
  cat("----- C-vine check (B) -----\n")
  
  cat("g =", g, "\n\n")
  
  
  
  cat("[1] vine$Matrix (structure):\n")
  
  print(vine$Matrix)
  
  
  
  cat("\n[2] summary(vine):\n")
  
  print(summary(vine))
  
  
  
  cat("\n[3] plot(vine): (Tree1 should look like star if C-vine)\n")
  
  invisible(try(plot(vine), silent = TRUE))
  
  
  
  invisible(vine)
  
}



cat("\nHow to check C-vine for B:\n",
    
    "Example: check_B_cvine(model_map_B, g=groups[1])\n")

