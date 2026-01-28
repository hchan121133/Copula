############################################################

# FINAL A) Increment-by-(g,t) "C-VINE" (3D) + Simulation

#

# - Each (g,t): fit hurdle margins for A_inc/B_inc/C_inc

# - Fit C-vine copula (forced) on PIT(U) (dimension=3)

# - Simulate increments -> cumulate -> exceed probability

# - Export:

#   (1) exceed summary CSV

#   (2) exceed-by-year CSV

#   (3) margin fit summary CSV (which dist used)

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

# [A-0] USER CONTROLS (여기만 바꾸면 운영 가능)

# =========================================================

LIMIT <- 10e8         # 10억 초과 여부 기준

Tmax  <- 20           # duration 1..20 가정



min_cell_n <- 300     # (g,t) 표본이 작으면 g풀링(fallback) 사용

min_pos_n  <- 80      # 양수부 표본이 작으면 ECDF로 대체(모수분포 적합 불안정 방지)



n_sim <- 1000        # ★ 시뮬레이션 횟수: 여기서 지정 ★



out_dir1 <- "result_Model1"      # 결과 저장 폴더



# =========================================================

# [A-1] INPUT CHECK

# =========================================================

stopifnot(exists("cum_df"))

need_cols <- c("id","g","duration","A_cum","B_cum","C_cum")

stopifnot(all(need_cols %in% names(cum_df)))



if (!dir.exists(out_dir1)) dir.create(out_dir1, recursive = TRUE)

stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")



# =========================================================

# [A-2] CUM -> INC (증분 생성: 단조성 보장 핵심)

# =========================================================

inc_df <- cum_df %>%
  
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

# [A-3] HURDLE MARGIN UTILS

# - p0 = P(X=0)

# - X|X>0 : lnorm/gamma/weibull/ECDF 중 AIC 최소 선택

# - randomized PIT 로 0 tie 문제 완화

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
  
  
  
  # 0 지급(또는 비정상)은 U~Unif(0,p0)로 랜덤화 → tie 완화
  
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

# [A-4] FIT CELL MODEL: "FORCE C-VINE" (dimension=3)

# =========================================================

fit_cell_model_cvine <- function(df_cell) {
  
  
  
  # (1) 마진 hurdle 적합
  
  mA <- fit_hurdle_margin(df_cell$A_inc, min_pos = min_pos_n)
  
  mB <- fit_hurdle_margin(df_cell$B_inc, min_pos = min_pos_n)
  
  mC <- fit_hurdle_margin(df_cell$C_inc, min_pos = min_pos_n)
  
  
  
  # (2) PIT → U(0,1)^3
  
  U <- cbind(
    
    pit_hurdle(df_cell$A_inc, mA),
    
    pit_hurdle(df_cell$B_inc, mB),
    
    pit_hurdle(df_cell$C_inc, mC)
    
  )
  
  
  
  # (3) family 후보 제한(속도/안정성)
  
  famset <- c(1, 3, 4, 5, 6, 13, 14, 16)
  
  
  
  # (4) ★ C-vine을 강제해서 적합 ★
  
  # - CVineStructureSelect는 "C-vine 구조" 전제 하에 pair-copula들을 선택
  
  vine <- VineCopula::RVineStructureSelect(
    
    U,
    
    familyset = famset,
    
    selectioncrit = "AIC",
    
    indeptest = TRUE,
    
    method = "mle",
    
    type = "CVine"
    
  )
  
  
  
  list(margins = list(A=mA, B=mB, C=mC), vine = vine)
  
}



# =========================================================

# [A-5] FIT ALL (g,t): 셀 부족 시 fallback

# =========================================================

groups <- sort(unique(inc_df$g))

Ts <- sort(unique(inc_df$duration))



handlers(global = TRUE)

plan(multisession, workers = 2)



model_map <- NULL



with_progress({
  
  p <- progressor(along = groups)
  
  
  
  model_map <- future_lapply(groups, function(gv) {
    
    df_g <- inc_df %>% dplyr::filter(g == gv)
    
    
    
    fallback <- fit_cell_model_cvine(df_g)
    
    
    
    by_t <- vector("list", length(Ts))
    
    names(by_t) <- as.character(Ts)
    
    
    
    for (t in Ts) {
      
      df_gt <- df_g %>% dplyr::filter(duration == t)
      
      if (nrow(df_gt) < min_cell_n) {
        
        by_t[[as.character(t)]] <- fallback
        
      } else {
        
        by_t[[as.character(t)]] <- fit_cell_model_cvine(df_gt)
        
      }
      
    }
    
    
    
    p(sprintf("fit A g=%s", gv))
    
    list(fallback = fallback, by_t = by_t)
    
  })
  
  
  
  names(model_map) <- groups
  
})



# =========================================================

# [A-6] SIMULATION

# =========================================================

simulate_group_paths_A <- function(gv, model_by_t, n_sim, Tmax, LIMIT, seed) {
  
  set.seed(seed)
  
  
  
  A_inc_mat <- matrix(0, nrow=n_sim, ncol=Tmax)
  
  B_inc_mat <- matrix(0, nrow=n_sim, ncol=Tmax)
  
  C_inc_mat <- matrix(0, nrow=n_sim, ncol=Tmax)
  
  
  
  for (t in 1:Tmax) {
    
    mdl <- model_by_t[[as.character(t)]]
    
    
    
    # C-vine에서 U 샘플링
    
    U <- RVineSim(n_sim, mdl$vine)
    
    
    
    # U → 금액(증분) 복원
    
    A_inc_mat[,t] <- pmax(q_hurdle(U[,1], mdl$margins$A), 0)
    
    B_inc_mat[,t] <- pmax(q_hurdle(U[,2], mdl$margins$B), 0)
    
    C_inc_mat[,t] <- pmax(q_hurdle(U[,3], mdl$margins$C), 0)
    
  }
  
  
  
  # 누적(단조성 자동 보장)
  
  A_cum <- t(apply(A_inc_mat, 1, cumsum))
  
  B_cum <- t(apply(B_inc_mat, 1, cumsum))
  
  C_cum <- t(apply(C_inc_mat, 1, cumsum))
  
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



with_progress({
  
  p <- progressor(along = groups)
  
  
  
  sim_res <- future_lapply(groups, function(gv) {
    
    out <- simulate_group_paths_A(
      
      gv = gv,
      
      model_by_t = model_map[[gv]]$by_t,
      
      n_sim = n_sim,
      
      Tmax = Tmax,
      
      LIMIT = LIMIT,
      
      seed = 1000 + match(gv, groups)
      
    )
    
    p(sprintf("sim A g=%s", gv))
    
    out
    
  })
  
})



# =========================================================

# [A-7] OUTPUT: EXCEED RESULTS

# =========================================================

A_exceed_summary <- dplyr::bind_rows(lapply(sim_res, function(x) {
  
  tibble(
    
    method = "A",
    
    g = x$g,
    
    n_sim = x$n_sim,
    
    exceed_prob_any = x$exceed_prob_any,
    
    S_cum_last_p50 = as.numeric(quantile(x$S_cum_last, 0.50)),
    
    S_cum_last_p90 = as.numeric(quantile(x$S_cum_last, 0.90)),
    
    S_cum_last_p99 = as.numeric(quantile(x$S_cum_last, 0.99))
    
  )
  
})) %>% dplyr::arrange(dplyr::desc(exceed_prob_any))



A_exceed_by_t <- dplyr::bind_rows(lapply(sim_res, function(x) {
  
  tibble(
    
    method = "A",
    
    g = x$g,
    
    duration = 1:Tmax,
    
    exceed_prob_t = x$exceed_prob_by_t
    
  )
  
}))



A_exceed_summary_path <- file.path(out_dir1, paste0("A_exceed_summary_", stamp, ".csv"))

A_exceed_by_t_path    <- file.path(out_dir1, paste0("A_exceed_by_t_", stamp, ".csv"))



write.csv(A_exceed_summary, A_exceed_summary_path, row.names = FALSE)

write.csv(A_exceed_by_t,    A_exceed_by_t_path,    row.names = FALSE)



cat("Saved(A) exceed:\n", A_exceed_summary_path, "\n", A_exceed_by_t_path, "\n")



# =========================================================

# [A-8] OUTPUT: "MARGIN FIT SUMMARY" CSV

# - 어떤 분포(lnorm/gamma/weibull/ecdf)로 잡혔는지 전체 덤프

# =========================================================

extract_A_margin_summary <- function(model_map) {
  
  all_g <- names(model_map)
  
  out_list <- list()
  
  idx <- 1
  
  
  
  for (g in all_g) {
    
    by_t <- model_map[[g]]$by_t
    
    all_t <- names(by_t)
    
    
    
    for (t in all_t) {
      
      mdl <- by_t[[t]]
      
      mlist <- mdl$margins
      
      
      
      for (k in names(mlist)) {
        
        m <- mlist[[k]]
        
        pos_type <- m$pos_model$type
        
        
        
        params <- NA_character_
        
        if (!is.null(m$pos_model$fit)) {
          
          est <- m$pos_model$fit$estimate
          
          params <- paste(names(est), round(as.numeric(est), 6), collapse="; ")
          
        }
        
        
        
        out_list[[idx]] <- data.frame(
          
          method = "A",
          
          g = g,
          
          t = as.integer(t),
          
          margin = k,       # A/B/C
          
          p0 = m$p0,
          
          pos_type = pos_type,
          
          params = params,
          
          stringsAsFactors = FALSE
          
        )
        
        idx <- idx + 1
        
      }
      
    }
    
  }
  
  dplyr::bind_rows(out_list) %>% dplyr::arrange(g, t, margin)
  
}



A_margin_tbl <- extract_A_margin_summary(model_map)

A_margin_path <- file.path(out_dir1, paste0("A_margin_fit_summary_", stamp, ".csv"))

write.csv(A_margin_tbl, A_margin_path, row.names = FALSE)

cat("Saved(A) margin fit summary:\n", A_margin_path, "\n")



# =========================================================

# [A-9] OUTPUT: "VINE FAMILY COUNTS" CSV

# - copula family(2변량 copula 타입) 빈도 요약

# =========================================================

family_name <- function(code){
  
  nm <- VineCopula::BiCopName(code)
  
  ifelse(is.na(nm), as.character(code), nm)
  
}



extract_A_family_counts <- function(model_map) {
  
  all_g <- names(model_map)
  
  out_list <- list()
  
  idx <- 1
  
  
  
  for (g in all_g) {
    
    by_t <- model_map[[g]]$by_t
    
    all_t <- names(by_t)
    
    
    
    for (t in all_t) {
      
      vine <- by_t[[t]]$vine
      
      fam_vec <- as.vector(vine$family)
      
      fam_vec <- fam_vec[fam_vec != 0]  # 0 제외(빈칸/대각/무의미)
      
      
      
      if (length(fam_vec) == 0) next
      
      
      
      tab <- sort(table(fam_vec), decreasing = TRUE)
      
      
      
      out_list[[idx]] <- data.frame(
        
        method = "A",
        
        g = g,
        
        t = as.integer(t),
        
        family_code = as.integer(names(tab)),
        
        family_name = family_name(as.integer(names(tab))),
        
        n_edges = as.integer(tab),
        
        stringsAsFactors = FALSE
        
      )
      
      idx <- idx + 1
      
    }
    
  }
  
  
  
  dplyr::bind_rows(out_list) %>%
    
    dplyr::arrange(g, t, dplyr::desc(n_edges))
  
}



A_family_tbl <- extract_A_family_counts(model_map)

A_family_path <- file.path(out_dir1, paste0("A_vine_family_counts_", stamp, ".csv"))

write.csv(A_family_tbl, A_family_path, row.names = FALSE)

cat("Saved(A) vine family counts:\n", A_family_path, "\n")



# =========================================================

# [A-10] "IS IT REALLY C-VINE?" 확인 방법

# =========================================================

# (1) 우리가 CVineStructureSelect로 적합했으므로 "구조 자체는 C-vine"이 전제입니다.

# (2) 그래도 확인하고 싶으면 아래를 보세요.

#

# - 특정 g,t의 vine 객체:

#     vine <- model_map[[g]]$by_t[[as.character(t)]]$vine

#

# - 확인 포인트:

#   a) print(vine$Matrix): C-vine은 1st tree에서 한 루트가 여러 변수와 star 형태로 연결

#   b) summary(vine): pair copula 목록/파라미터 출력

#   c) plot(vine): tree 1 구조가 star면 C-vine 직관 확인 가능



check_A_cvine <- function(model_map, g, t) {
  
  vine <- model_map[[g]]$by_t[[as.character(t)]]$vine
  
  cat("----- C-vine check (A) -----\n")
  
  cat("g =", g, "t =", t, "\n\n")
  
  
  
  cat("[1] vine$Matrix (structure):\n")
  
  print(vine$Matrix)
  
  
  
  cat("\n[2] summary(vine):\n")
  
  print(summary(vine))
  
  
  
  cat("\n[3] plot(vine): (Tree1 should look like star if C-vine)\n")
  
  # plot은 그래픽 디바이스 필요
  
  invisible(try(plot(vine), silent = TRUE))
  
  
  
  invisible(vine)
  
}



cat("\nHow to check C-vine for A:\n",
    
    "Example: check_A_cvine(model_map, g=groups[1], t=1)\n")







