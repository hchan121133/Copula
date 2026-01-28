## =============================================================================
## 04_evaluation.R : 결과물 평가 (A vs B vs C)
## =============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2)
})

cat(">>> [4/5] 모델 결과 비교 (Evaluation)...\n")

res_list <- list()
if(file.exists("output/res_method_A.rds")) res_list[["A"]] <- readRDS("output/res_method_A.rds")
if(file.exists("output/res_method_B.rds")) res_list[["B"]] <- readRDS("output/res_method_B.rds")
if(file.exists("output/res_method_C.rds")) res_list[["C"]] <- readRDS("output/res_method_C.rds")

if(length(res_list) > 0) {
  all_res <- bind_rows(res_list)
  
  # 1. 20차년도 최종 탈퇴율(한도초과율) 비교 테이블
  comp_table <- all_res %>%
    filter(Duration == 20) %>%
    select(Method, Group, Exceed_Prob) %>%
    pivot_wider(names_from=Method, values_from=Exceed_Prob)
  
  # 차이 계산 (B가 있으면 기준)
  if("B" %in% names(comp_table) && "C" %in% names(comp_table)) {
    comp_table$Diff_C_B <- comp_table$C - comp_table$B
  }
  
  write.csv(comp_table, "output/04_Model_Comparison_ABC.csv", row.names=FALSE)
  print(head(comp_table))
  
  # 2. 시각화 (Duration별 변화 추세 - 첫번째 그룹만 예시)
  target_g <- comp_table$Group[1]
  plot_data <- all_res %>% filter(Group == target_g)
  
  p <- ggplot(plot_data, aes(x=Duration, y=Exceed_Prob, color=Method)) +
    geom_line(size=1.2) +
    geom_point() +
    theme_minimal() +
    labs(title=paste("Exceed Probability by Method (Group:", target_g, ")"),
         y="Prob(Sum > 10e8)")
  
  ggsave("output/04_Method_Comparison_Plot.png", p, width=8, height=6)
  
  cat(">>> 완료: output/04_Model_Comparison_ABC.csv 및 Plot 저장됨\n")
} else {
  cat(">>> 경고: 비교할 결과 파일이 없습니다.\n")
}