

############################################################

# cvine_cum.R  (METHOD = CUM)  --- DATA GEN FINAL (No CAT)

#

# 목표:

# - "현재 담보( CC / SH / HP / SG ) 안에서" 상관관계(동시·동연도 고액)를 강화

# - 20년 누적 10억(LIMIT) 초과 케이스를 유의미하게 증가

# - raw_df / cum_df 생성 후 CSV(.csv.gz)로 Export

#

# 핵심 튜닝(상관강화):

# 1) severe(공통 중증상태) + Z_persist(지속성) 추가

# 2) CC 고액 트리거(CC3/CC5) → SH/HP/SG를 조건부로 끌어올림

# 3) HP cap 완화(30→60) : 단가 동일, "일수(빈도)"로 누적가속

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

# [DATA GEN] 가상 데이터 생성 (상관강화 + 10억 초과 튜닝)

# a = CC1~CC8 + SH1~SH4, b = HP, c = SG1~SG5

############################################################

if (!exists("raw_df") || !exists("cum_df")) {
  
  
  
  set.seed(20260206)
  
  
  
  # -------------------------
  
  # CONTROL / TUNING KNOBS
  
  # -------------------------
  
  N     <- 20000
  
  Tmax  <- 20
  
  LIMIT <- 10e8
  
  
  
  # severe(공통충격) 강도/발생률
  
  SEV_INTERCEPT <- -2.6  # (더 크게: -2.4, -2.2 ...) => 중증↑, 10억초과↑
  
  SEV_AGE_SLOPE <-  0.06
  
  SEV_LOADING   <-  1.45 # 공통요인 로딩↑ => 동시발생↑
  
  PERSIST_W     <-  0.75 # 지속성↑ => 연속클러스터↑
  
  
  
  # 트리거(고액) 연동 강도
  
  TRG_CC5_to_SH <-  0.90
  
  TRG_CC5_to_HP <-  0.85
  
  TRG_CC5_to_SG <-  0.90
  
  
  
  # HP 일수 cap (단가 동일, 빈도 누적으로 가속)
  
  HP_CAP_DAYS   <- 60
  
  
  
  # -------------------------
  
  # AMOUNTS
  
  # -------------------------
  
  amt_CC1 <- 1000  * 10^4
  
  amt_CC2 <- 200   * 10^4
  
  amt_CC3 <- 5000  * 10^4
  
  amt_CC4 <- 1000  * 10^4
  
  amt_CC5 <- 10000 * 10^4
  
  amt_CC6 <- 2000  * 10^4
  
  amt_CC7 <- 1000  * 10^4
  
  amt_CC8 <- 200   * 10^4
  
  
  
  amt_SH1 <- 2000 * 10^4
  
  amt_SH2 <- 2000 * 10^4
  
  amt_SH3 <- 2000 * 10^4
  
  amt_SH4 <- 2000 * 10^4
  
  
  
  amt_HP <- 10 * 10^4
  
  
  
  amt_SG <- c(10 * 10^4, 15 * 10^4, 100 * 10^4, 250 * 10^4, 500 * 10^4)
  
  
  
  CAP_Sick_SG_PER_YEAR <- 8
  
  
  
  logit <- function(x) 1/(1+exp(-x))
  
  
  
  # -------------------------
  
  # ID / LATENT
  
  # -------------------------
  
  id_tbl <- tibble(
    
    id        = sprintf("ID%06d", 1:N),
    
    sex       = sample(c("M","F"), N, TRUE),
    
    issue_age = sample(15:75, N, TRUE)
    
  ) %>%
    
    mutate(
      
      age_band = cut(
        
        issue_age,
        
        breaks = c(0,19,29,39,49,59,69,120),
        
        labels = c("15-19","20-29","30-39","40-49","50-59","60-69","70+")
        
      )
      
    )
  
  
  
  latent_tbl <- tibble(
    
    id       = id_tbl$id,
    
    Z_common = rnorm(N, sd = 1.7),  # tail 쪽 상관 강화
    
    Z_HP     = rnorm(N, sd = 1.3),
    
    Z_surg   = rnorm(N, sd = 1.3)
    
  )
  
  
  
  panel <- id_tbl %>%
    
    crossing(duration = 1:Tmax) %>%
    
    mutate(attained_age = issue_age + duration - 1) %>%
    
    left_join(latent_tbl, by = "id") %>%
    
    group_by(id) %>%
    
    mutate(
      
      # 지속성(연도별 흔들림을 조금 주되 id별로 클러스터링)
      
      Z_persist = PERSIST_W*Z_common + (1-PERSIST_W)*rnorm(n()),
      
      
      
      # 공통 중증상태(severe): 동시발생/동시고액을 만드는 이산 충격
      
      severe_p = logit(SEV_INTERCEPT + SEV_AGE_SLOPE*(attained_age-50) + SEV_LOADING*Z_persist),
      
      severe   = rbinom(n(), 1, prob = severe_p)
      
    ) %>%
    
    ungroup()
  
  
  
  # -------------------------
  
  # CC (고액 트리거 포함)
  
  # -------------------------
  
  cc <- panel %>%
    
    mutate(
      
      # 빈도성(포아송): severe면 log-lambda가 올라가 동시다발↑
      
      lambda_cc1 = exp(-4.6 + 0.05*(attained_age-40) + 0.03*(duration-1) + 0.85*Z_common + 0.60*severe),
      
      CC_n1 = pmin(rpois(n(), lambda = lambda_cc1), CAP_Sick_SG_PER_YEAR),
      
      CC_y1 = CC_n1 * amt_CC1,
      
      
      
      lambda_cc2 = exp(-4.6 + 0.05*(attained_age-40) + 0.03*(duration-1) + 0.75*Z_common + 0.55*severe),
      
      CC_n2 = pmin(rpois(n(), lambda = lambda_cc2), CAP_Sick_SG_PER_YEAR),
      
      CC_y2 = CC_n2 * amt_CC2,
      
      
      
      # 고액 트리거(이항): severe 강하게 반영
      
      CC_n3 = rbinom(n(), 1, prob = logit(-6.0 + 0.06*(attained_age-40) + 0.80*Z_common + 1.15*severe)),
      
      CC_y3 = CC_n3 * amt_CC3,
      
      
      
      CC_n4 = rbinom(n(), 1, prob = logit(-5.8 + 0.06*(attained_age-40) + 0.75*Z_common + 0.85*severe)),
      
      CC_y4 = CC_n4 * amt_CC4,
      
      
      
      CC_n5 = rbinom(n(), 1, prob = logit(-5.6 + 0.07*(attained_age-40) + 0.85*Z_common + 1.35*severe)),
      
      CC_y5 = CC_n5 * amt_CC5,
      
      
      
      CC_n6 = rbinom(n(), 1, prob = logit(-5.9 + 0.07*(attained_age-45) + 0.75*Z_common + 0.75*severe)),
      
      CC_y6 = CC_n6 * amt_CC6,
      
      
      
      CC_n7 = rbinom(n(), 1, prob = logit(-5.6 + 0.08*(attained_age-45) + 0.30*(sex=="M") + 0.85*Z_common + 0.80*severe)),
      
      CC_y7 = CC_n7 * amt_CC7,
      
      
      
      CC_n8 = rbinom(n(), 1, prob = logit(-5.9 + 0.08*(attained_age-45) + 0.25*(sex=="M") + 0.90*Z_common + 0.85*severe)),
      
      CC_y8 = CC_n8 * amt_CC8,
      
      
      
      CC_pay   = CC_y1 + CC_y2 + CC_y3 + CC_y4 + CC_y5 + CC_y6 + CC_y7 + CC_y8,
      
      CC_event = as.integer(CC_pay > 0),
      
      
      
      CC3_event = as.integer(CC_n3 > 0),
      
      CC5_event = as.integer(CC_n5 > 0)
      
    ) %>%
    
    dplyr::select(id, duration,CC_n1:CC_n8, CC_y1:CC_y8, CC_pay, CC_event, CC3_event, CC5_event)
  
  
  
  # -------------------------
  
  # SH (CC5 트리거 연동)
  
  # -------------------------
  
  sh <- panel %>%
    
    left_join(cc %>% dplyr::select(id, duration, CC3_event, CC5_event), by = c("id","duration")) %>%
    
    mutate(
      
      lambda_sh1 = exp(-4.6 + 0.05*(attained_age-40) + 0.03*(duration-1) + 0.85*Z_common + 0.55*severe + 0.25*CC5_event),
      
      SH_n1 = pmin(rpois(n(), lambda = lambda_sh1), CAP_Sick_SG_PER_YEAR),
      
      SH_y1 = SH_n1 * amt_SH1,
      
      
      
      lambda_sh2 = exp(-4.6 + 0.05*(attained_age-40) + 0.03*(duration-1) + 0.75*Z_common + 0.50*severe + 0.20*CC5_event),
      
      SH_n2 = pmin(rpois(n(), lambda = lambda_sh2), CAP_Sick_SG_PER_YEAR),
      
      SH_y2 = SH_n2 * amt_SH2,
      
      
      
      SH_n3 = rbinom(n(), 1, prob = logit(-5.6 + 0.08*(attained_age-45) + 0.85*Z_common + 0.95*severe + TRG_CC5_to_SH*CC5_event)),
      
      SH_y3 = SH_n3 * amt_SH3,
      
      
      
      SH_n4 = rbinom(n(), 1, prob = logit(-5.9 + 0.08*(attained_age-45) + 0.90*Z_common + 1.00*severe + (TRG_CC5_to_SH+0.10)*CC5_event)),
      
      SH_y4 = SH_n4 * amt_SH4,
      
      
      
      SH_pay   = SH_y1 + SH_y2 + SH_y3 + SH_y4,
      
      SH_event = as.integer(SH_pay > 0)
      
    ) %>%
    
    dplyr::select(id, duration,
                  
                  SH_n1:SH_n4, SH_y1:SH_y4, SH_pay, SH_event)
  
  
  
  # -------------------------
  
  # Combine + HP
  
  # -------------------------
  
  tmp <- panel %>%
    
    left_join(cc, by = c("id","duration")) %>%
    
    left_join(sh, by = c("id","duration")) %>%
    
    mutate(
      
      # HP admit 확률: severe + CC/SH + CC5 트리거로 강하게 동조
      
      HP_admit_p = logit(
        
        -3.4 + 0.04*(attained_age-40) +
          
          0.65*Z_common +
          
          1.60*CC_event + 1.60*SH_event +
          
          1.05*severe + TRG_CC5_to_HP*CC5_event
        
      ),
      
      HP_admit = rbinom(n(), 1, prob = HP_admit_p),
      
      
      
      # 입원일수(=빈도): severe/CC5이면 크게 증가
      
      HP_days_mu = pmax(1, 6 + 6*CC_event + 6*SH_event + 2.0*Z_common + 10*severe + 7*CC5_event),
      
      HP_n = ifelse(HP_admit==1, rpois(n(), lambda = HP_days_mu), 0),
      
      HP_n = pmin(HP_n, HP_CAP_DAYS),
      
      
      
      HP_pay   = HP_n * amt_HP,
      
      HP_event = as.integer(HP_pay > 0)
      
    )
  
  
  
  # -------------------------
  
  # SG (수술) : high-class(4,5)에 severe/CC5/HP 연동 강화
  
  # -------------------------
  
  for (k in 1:5) {
    
    base_p <- 0.020 + 0.003*k
    
    high_class_boost <- ifelse(k >= 4, 0.95, 0.00)  # (기존 0.60 수준보다 강하게)
    
    
    
    p_k <- logit(
      
      qlogis(base_p) +
        
        0.70*tmp$Z_surg +
        
        0.55*tmp$Z_common +
        
        1.40*tmp$CC_event +
        
        1.40*tmp$HP_event +
        
        0.95*tmp$severe +
        
        TRG_CC5_to_SG*tmp$CC5_event +
        
        high_class_boost
      
    )
    
    
    
    tmp[[paste0("SG_n", k)]] <- rbinom(nrow(tmp), 1, prob = p_k)
    
    tmp[[paste0("SG_y", k)]] <- tmp[[paste0("SG_n", k)]] * amt_SG[k]
    
  }
  
  
  
  tmp <- tmp %>%
    
    mutate(
      
      SG_pay = SG_y1 + SG_y2 + SG_y3 + SG_y4 + SG_y5,
      
      pay_total = CC_pay + SH_pay + HP_pay + SG_pay
      
    )
  
  
  
  # -------------------------
  
  # raw_df
  
  # -------------------------
  
  raw_df <- tmp %>%
    
    arrange(id, duration) %>%
    
    group_by(id) %>%
    
    mutate(
      
      S_cum       = cumsum(pay_total),
      
      exceed_10e8 = (S_cum > LIMIT)
      
    ) %>%
    
    ungroup() %>%
    
    dplyr::select(
      
      id, sex, issue_age, age_band, duration, attained_age,
      
      severe, Z_common, Z_persist,
      
      CC_n1:CC_n8, CC_y1:CC_y8, CC_pay, CC_event, CC3_event, CC5_event,
      
      SH_n1:SH_n4, SH_y1:SH_y4, SH_pay, SH_event,
      
      HP_n, HP_pay, HP_event,
      
      SG_n1:SG_n5, SG_y1:SG_y5, SG_pay,
      
      pay_total, S_cum, exceed_10e8
      
    )
  
  
  
  # -------------------------
  
  # quick diagnostics
  
  # -------------------------
  
  exceed_id_rate <- raw_df %>%
    
    group_by(id) %>%
    
    summarise(exceed_any = any(S_cum > LIMIT), .groups = "drop") %>%
    
    summarise(rate = mean(exceed_any), n_exceed = sum(exceed_any))
  
  print(exceed_id_rate)
  
  
  
  check_tbl <- raw_df %>%
    
    mutate(main_event = (CC_pay + SH_pay > 0)) %>%
    
    summarise(
      
      HP_rate_when_a     = mean(HP_pay > 0 & main_event),
      
      HP_rate_when_not_a = mean(HP_pay > 0 & !main_event),
      
      SG_rate_when_a     = mean(SG_pay > 0 & main_event),
      
      SG_rate_when_not_a = mean(SG_pay > 0 & !main_event),
      
      severe_rate        = mean(severe==1)
      
    )
  
  print(check_tbl)
  
  
  
  # -------------------------
  
  # cum_df (CUM Method용)
  
  # -------------------------
  
  cum_df <- raw_df %>%
    
    arrange(id, duration) %>%
    
    group_by(id) %>%
    
    mutate(
      
      A_cum = cumsum(CC_pay) + cumsum(SH_pay),
      
      B_cum = cumsum(HP_pay),
      
      C_cum = cumsum(SG_pay),
      
      S_cum_abc = A_cum + B_cum + C_cum
      
    ) %>%
    
    ungroup() %>%
    
    mutate(g = paste0(sex, "_", age_band))
  
  
  
  ############################################################
  
  # [EXPORT] CSV(.gz)로 저장
  
  ############################################################
  
  out_dir <- "result_RAW"
  
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  
  
  stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  
  raw_path <- file.path(out_dir, paste0("raw_df_", stamp, ".csv"))
  
  cum_path <- file.path(out_dir, paste0("cum_df_", stamp, ".csv"))
  
  
  
  write.csv(raw_df, raw_path, row.names = FALSE)
  
  write.csv(cum_df, cum_path, row.names = FALSE)
  
  
  
  cat("Saved:\n", raw_path, "\n", cum_path, "\n")
  
}



