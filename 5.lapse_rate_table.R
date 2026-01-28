## =============================================================================
## 05_lapse_rate_table.R : 최종 탈퇴율 테이블 생성 (Method 선택 가능)
## =============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr)
})

# ====================================================
# [User Option] 사용할 모델 선택 ("A", "B", "C")
SELECTED_METHOD <- "C"  # <-- 여기서 변경하세요
# ====================================================

file_map <- c("A"="output/res_method_A.rds", 
              "B"="output/res_method_B.rds", 
              "C"="output/res_method_C.rds")

target_file <- file_map[SELECTED_METHOD]

if(!file.exists(target_file)) {
  stop(paste("선택한 Method", SELECTED_METHOD, "의 결과 파일이 없습니다."))
}

cat(">>> [5/5] 최종 산출물 생성 (Selected Method:", SELECTED_METHOD, ")\n")
res_df <- readRDS(target_file)

# g 컬럼 분리 및 Pivot
final_wide <- res_df %>%
  separate(Group, into=c("Sex", "Age_Band"), sep="_") %>%
  select(Sex, Age_Band, Duration, Exceed_Prob) %>%
  arrange(Sex, Age_Band, Duration) %>%
  pivot_wider(
    names_from = Duration,
    values_from = Exceed_Prob,
    names_prefix = "Year_"
  )

out_name <- paste0("output/Final_Lapse_Rates_Method_", SELECTED_METHOD, ".csv")
write.csv(final_wide, out_name, row.names=FALSE)

cat("=======================================================\n")
cat(" [SUCCESS] 최종 산출 완료\n")
cat(" 저장 파일:", out_name, "\n")
cat("=======================================================\n")