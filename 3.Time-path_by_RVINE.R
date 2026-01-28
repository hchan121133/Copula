############################################################
# FINAL C) Time-path (3*T=60D) "R-VINE" + Simulation
#
# - For each group g:
#   1) Build 60D vector per id: A_inc_1..T, B_inc_1..T, C_inc_1..T
#   2) Fit hurdle margins for each dimension
#   3) Fit R-vine copula on 60D U (Structure Selection)
#   4) Simulate 60D paths -> cumulate -> exceed probability
# - Export:
#   (1) exceed summary CSV (C_exceed...)
#   (2) exceed-by-year CSV (C_exceed...)
#   (3) margin fit summary CSV
#   (4) vine family counts CSV
# - Includes "How to check structure" section
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
# [C-0] USER CONTROLS
# =========================================================
LIMIT <- 10e8
Tmax  <- 20
min_pos_n <- 80

n_sim <- 1000          # ★ 시뮬레이션 횟수 (안정성 위해 1000회 이상 권장)

TRUNC_LEVEL <- 4       # ★ 4단계까지 확장 (A1->A2 연결 강화)
# 고액 리스크 집중 Family Set
familyset   <- c(1, 3, 4, 5, 6, 13, 14, 16, 23, 24, 26, 33, 34, 36)

# ★ [수정] 평가 코드(evaluation.R)가 인식할 수 있도록 폴더명 변경
out_dir3 <- "result_Model3"

# =========================================================
# [C-1] INPUT CHECK
# =========================================================
stopifnot(exists("cum_df"))
need_cols <- c("id","g","duration","A_cum","B_cum","C_cum")
stopifnot(all(need_cols %in% names(cum_df)))

if (!dir.exists(out_dir3)) dir.create(out_dir3, recursive = TRUE)
stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

# =========================================================
# [C-2] CUM -> INC (long)
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
# [C-3] HURDLE UTILS (A와 동일)
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
# [C-4] BUILD 60D WIDE TABLE PER GROUP
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
# [C-5] FIT 60D "R-VINE" PER GROUP
# =========================================================
fit_timevine_group_rvine <- function(df_wide, min_pos_n, familyset, trunc_level) {
  
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
  
  # (3) R-vine 자동 구조 탐색 및 적합
  vine <- NULL
  
  vine <- tryCatch({
    RVineStructureSelect(
      data = U,
      familyset = familyset,
      type = 0,              # 0: R-Vine (General)
      selectioncrit = "AIC",
      indeptest = TRUE,
      level = 0.05,
      trunclevel = trunc_level, 
      progress = FALSE,
      method = "mle"
    )
  }, error = function(e) {
    message("Vine fit error (will skip this group): ", e$message)
    NULL
  })
  
  list(margins = margins, vine = vine, colnames = colnames(U))
}

# =========================================================
# [C-6] SIMULATE 60D PATHS
# =========================================================
simulate_timevine_group_C <- function(gv, mdl, n_sim, Tmax, LIMIT, seed) {
  set.seed(seed)
  
  if(is.null(mdl$vine)) return(NULL)
  
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
# [C-7] RUN FIT + SIM PER GROUP
# =========================================================
groups <- sort(unique(inc_long$g))

handlers(global = TRUE)
plan(multisession, workers = 2)

model_map_C <- NULL
sim_res_C <- NULL

with_progress({
  p <- progressor(along = groups)
  
  model_map_C <- future_lapply(groups, function(gv) {
    df_g <- inc_long %>% dplyr::filter(g == gv)
    df_wide <- make_wide_inc(df_g, Tmax = Tmax)
    
    mdl <- fit_timevine_group_rvine(
      df_wide = df_wide,
      min_pos_n = min_pos_n,
      familyset = familyset,
      trunc_level = TRUNC_LEVEL
    )
    
    p(sprintf("fit C g=%s", gv))
    mdl
  })
  
  names(model_map_C) <- groups
})

with_progress({
  p <- progressor(along = groups)
  
  sim_res_C <- future_lapply(groups, function(gv) {
    if(is.null(model_map_C[[gv]]$vine)) return(NULL)
    
    out <- simulate_timevine_group_C(
      gv = gv,
      mdl = model_map_C[[gv]],
      n_sim = n_sim,
      Tmax = Tmax,
      LIMIT = LIMIT,
      seed = 5000 + match(gv, groups)
    )
    p(sprintf("sim C g=%s", gv))
    out
  })
})

# NULL(실패한 그룹) 제거
sim_res_C <- sim_res_C[!sapply(sim_res_C, is.null)]

# 유효한 결과가 하나도 없는 경우 중단
if(length(sim_res_C) == 0) {
  stop("모든 그룹에 대한 시뮬레이션이 실패했습니다. 데이터나 파라미터를 확인해주세요.")
}

# =========================================================
# [C-8] OUTPUT: EXCEED RESULTS
# =========================================================
C_exceed_summary <- dplyr::bind_rows(lapply(sim_res_C, function(x) {
  tibble(
    method = "C",  # ★ [수정] Method 이름을 C로 설정
    g = x$g,
    n_sim = x$n_sim,
    exceed_prob_any = x$exceed_prob_any,
    S_cum_last_p50 = as.numeric(quantile(x$S_cum_last, 0.50)),
    S_cum_last_p90 = as.numeric(quantile(x$S_cum_last, 0.90)),
    S_cum_last_p99 = as.numeric(quantile(x$S_cum_last, 0.99))
  )
})) %>% dplyr::arrange(dplyr::desc(exceed_prob_any))

C_exceed_by_t <- dplyr::bind_rows(lapply(sim_res_C, function(x) {
  tibble(
    method = "C",  # ★ [수정] Method 이름을 C로 설정
    g = x$g,
    duration = 1:Tmax,
    exceed_prob_t = x$exceed_prob_by_t
  )
}))

# ★ [수정] 파일명 C_exceed_... 로 변경
C_exceed_summary_path <- file.path(out_dir3, paste0("C_exceed_summary_", stamp, ".csv"))
C_exceed_by_t_path    <- file.path(out_dir3, paste0("C_exceed_by_t_", stamp, ".csv"))

write.csv(C_exceed_summary, C_exceed_summary_path, row.names = FALSE)
write.csv(C_exceed_by_t,    C_exceed_by_t_path,    row.names = FALSE)

cat("Saved(C) exceed:\n", C_exceed_summary_path, "\n", C_exceed_by_t_path, "\n")

# =========================================================
# [C-9] OUTPUT: "MARGIN FIT SUMMARY" CSV (60D)
# =========================================================
extract_C_margin_summary <- function(model_map_C) {
  all_g <- names(model_map_C)
  out_list <- list()
  idx <- 1
  
  for (g in all_g) {
    mdl <- model_map_C[[g]]
    if(is.null(mdl$vine)) next
    
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
        method = "C", # ★ [수정]
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

C_margin_tbl <- extract_C_margin_summary(model_map_C)
C_margin_path <- file.path(out_dir3, paste0("C_margin_fit_summary_", stamp, ".csv"))
write.csv(C_margin_tbl, C_margin_path, row.names = FALSE)
cat("Saved(C) margin fit summary:\n", C_margin_path, "\n")

# =========================================================
# [C-10] OUTPUT: "VINE FAMILY COUNTS" CSV
# =========================================================
family_name <- function(code){
  nm <- VineCopula::BiCopName(code)
  ifelse(is.na(nm), as.character(code), nm)
}

extract_C_family_counts <- function(model_map_C) {
  all_g <- names(model_map_C)
  out_list <- list()
  idx <- 1
  
  for (g in all_g) {
    if(is.null(model_map_C[[g]]$vine)) next
    vine <- model_map_C[[g]]$vine
    
    fam_vec <- as.vector(vine$family)
    fam_vec <- fam_vec[fam_vec != 0]
    if (length(fam_vec) == 0) next
    
    tab <- sort(table(fam_vec), decreasing = TRUE)
    
    out_list[[idx]] <- data.frame(
      method = "C", # ★ [수정]
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

C_family_tbl <- extract_C_family_counts(model_map_C)
C_family_path <- file.path(out_dir3, paste0("C_vine_family_counts_", stamp, ".csv"))
write.csv(C_family_tbl, C_family_path, row.names = FALSE)
cat("Saved(C) vine family counts:\n", C_family_path, "\n")

# =========================================================
# [C-11] "IS IT REALLY R-VINE?" 확인 방법
# =========================================================
check_C_rvine <- function(model_map_C, g) {
  if(is.null(model_map_C[[g]]$vine)) {
    cat("Model for group", g, "failed to fit.\n")
    return(NULL)
  }
  
  vine <- model_map_C[[g]]$vine
  cat("----- R-vine check (C) -----\n")
  cat("g =", g, "\n\n")
  
  cat("[1] vine$Matrix (Structure):\n")
  print(vine$Matrix)
  
  cat("\n[2] summary(vine):\n")
  print(summary(vine))
  
  cat("\n[3] plot(vine): (Should show optimized tree structure)\n")
  invisible(try(plot(vine), silent = TRUE))
  
  invisible(vine)
}

cat("\nHow to check R-vine for C:\n",
    "Example: check_C_rvine(model_map_C, g=groups[1])\n")
