## =============================================================================
## 04_evaluation.R : 결과물 평가 (CSV 기반 A vs B 자동 비교)
## =============================================================================
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
})

cat(">>> [4/5] 모델 결과 비교 (Evaluation - CSV Mode)...\n")

# ------------------------------------------------------------------------------
# 1. 가장 최근 CSV 파일을 찾는 함수
# ------------------------------------------------------------------------------
get_latest_result <- function(folder, pattern, method_label) {
  if (!dir.exists(folder)) return(NULL)
  
  files <- list.files(folder, pattern = pattern, full.names = TRUE)
  if (length(files) == 0) return(NULL)
  
  # 가장 최근 파일 선택
  latest_file <- files[which.max(file.info(files)$mtime)]
  cat(sprintf("   -> Load [%s]: %s\n", method_label, basename(latest_file)))
  
  df <- read.csv(latest_file)
  
  # 컬럼명 표준화 (코드 호환성 유지)
  # CSV: method, g, duration, exceed_prob_t
  # Target: Method, Group, Duration, Exceed_Prob
  df %>%
    dplyr::mutate(Method = method_label) %>%
    dplyr::rename(
      Group = g,
      Duration = duration,
      Exceed_Prob = exceed_prob_t
    ) %>%
    dplyr::select(Method, Group, Duration, Exceed_Prob)
}

# ------------------------------------------------------------------------------
# 2. 결과 로드 (Model A, B, C 등)
# ------------------------------------------------------------------------------
res_list <- list()

# Model A (result_Model1 폴더의 A_exceed_by_t_*.csv)
res_list[["A"]] <- get_latest_result("result_Model1", "A_exceed_by_t_.*\\.csv", "A")

# Model B (result_Model3 폴더의 B_exceed_by_t_*.csv -> R-Vine 결과)
# ※ 만약 Model2(C-Vine)를 썼다면 result_Model2로 변경
res_list[["B"]] <- get_latest_result("result_Model2", "B_exceed_by_t_.*\\.csv", "B")

# Model C (옵션)
res_list[["C"]] <- get_latest_result("result_Model3", "C_exceed_by_t_.*\\.csv", "C")

# NULL 제거
res_list <- res_list[!sapply(res_list, is.null)]

# ------------------------------------------------------------------------------
# 3. 비교 및 시각화
# ------------------------------------------------------------------------------
if (length(res_list) > 0) {
  all_res <- dplyr::bind_rows(res_list)
  
  # (1) 저장 폴더 확인
  if (!dir.exists("evaluation")) dir.create("evaluation")
  
  # (2) 20차년도 최종 탈퇴율(한도초과율) 비교 테이블
  comp_table <- all_res %>%
    dplyr::filter(Duration == 20) %>%
    dplyr::select(Method, Group, Exceed_Prob) %>%
    tidyr::pivot_wider(names_from = Method, values_from = Exceed_Prob)
  
  # 차이 계산 (B가 기준)
  if ("A" %in% names(comp_table) && "B" %in% names(comp_table)) {
    comp_table$Diff_B_A <- comp_table$B - comp_table$A
  }
  
  out_csv <- "evaluation/04_Model_Comparison_ABC.csv"
  write.csv(comp_table, out_csv, row.names = FALSE)
  
  cat("\n[Comparison Table (Top 5)]\n")
  print(head(comp_table))
  
  # (3) 시각화 (Duration별 변화 추세 - 첫번째 그룹 예시)
  if (nrow(comp_table) > 0) {
    target_g <- comp_table$Group[1]
    plot_data <- all_res %>% dplyr::filter(Group == target_g)
    
    p <- ggplot(plot_data, aes(x = Duration, y = Exceed_Prob, color = Method)) +
      geom_line(linewidth = 1.2) +
      geom_point(size = 2) +
      theme_minimal() +
      labs(
        title = paste("Exceed Probability Path (Group:", target_g, ")"),
        subtitle = "Model A (Yearly) vs Model B (Time-Path/R-Vine)",
        y = "Prob(Sum > 10e8)",
        x = "Duration (Year)"
      )
    
    out_plot <- "evaluation/04_Method_Comparison_Plot.png"
    ggsave(out_plot, p, width = 8, height = 6)
    cat(sprintf("\n>>> 완료: \n 1. %s\n 2. %s\n", out_csv, out_plot))
    
  } else {
    cat(">>> 데이터가 없어 그래프를 그릴 수 없습니다.\n")
  }
  
} else {
  cat(">>> [경고] 비교할 결과 파일(.csv)을 찾지 못했습니다.\n")
  cat("    result_Model1 또는 result_Model3 폴더에 csv 파일이 있는지 확인해주세요.\n")
}
