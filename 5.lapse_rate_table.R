## =============================================================================
## 05_lapse_rate_table.R : 최종 탈퇴율 테이블 생성 (CSV 자동 감지)
## =============================================================================
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
})

# ====================================================
# [User Option] 사용할 모델 선택 ("A", "B", "C")
# A: result_Model1 폴더 (Increment)
# B: result_Model3 폴더 (R-Vine)
# C: result_ModelC 폴더 (Cumulative R-Vine)
# ====================================================
SELECTED_METHOD <- "B"  # <-- 사용할 모델 (A, B, C 중 택 1)

# ------------------------------------------------------------------------------
# 1. 설정: 모델별 저장 폴더 및 파일 패턴 매핑
# ------------------------------------------------------------------------------
config_map <- list(
  "A" = list(dir = "result_Model1", pattern = "A_exceed_by_t_.*\\.csv"),
  "B" = list(dir = "result_Model2", pattern = "B_exceed_by_t_.*\\.csv"),
  "C" = list(dir = "result_Model3", pattern = "C_exceed_by_t_.*\\.csv")
)

target_cfg <- config_map[[SELECTED_METHOD]]
if (is.null(target_cfg)) stop("잘못된 Method 선택입니다. A, B, C 중 하나를 선택하세요.")

# ------------------------------------------------------------------------------
# 2. 최신 결과 파일 찾기
# ------------------------------------------------------------------------------
cat(sprintf(">>> [5/5] 최종 산출물 생성 (Method: %s)\n", SELECTED_METHOD))
cat(sprintf("    검색 폴더: %s\n", target_cfg$dir))

if (!dir.exists(target_cfg$dir)) {
  stop(sprintf("폴더 '%s'가 존재하지 않습니다. 해당 모델 시뮬레이션을 먼저 실행하세요.", target_cfg$dir))
}

files <- list.files(target_cfg$dir, pattern = target_cfg$pattern, full.names = TRUE)
if (length(files) == 0) {
  stop(sprintf("폴더 내에 패턴 '%s'에 맞는 결과 파일이 없습니다.", target_cfg$pattern))
}

# 가장 최근 수정된 파일 선택
latest_file <- files[which.max(file.info(files)$mtime)]
cat(sprintf("    사용 파일: %s\n", basename(latest_file)))

# ------------------------------------------------------------------------------
# 3. 데이터 로드 및 변환
# ------------------------------------------------------------------------------
res_df <- read.csv(latest_file)

# 컬럼명 확인 및 표준화 (소문자 g, duration 등을 대문자로 변환)
# CSV 표준 출력: method, g, duration, exceed_prob_t
if (!"g" %in% names(res_df) | !"duration" %in% names(res_df)) {
  stop("CSV 파일의 컬럼명이 예상과 다릅니다. (g, duration 컬럼 필요)")
}

final_wide <- res_df %>%
  # 필요한 컬럼만 선택 및 이름 변경
  dplyr::select(Group = g, Duration = duration, Exceed_Prob = exceed_prob_t) %>%
  # Group(예: M_30-39)을 성별/연령대로 분리
  tidyr::separate(Group, into = c("Sex", "Age_Band"), sep = "_") %>%
  # 정렬
  dplyr::arrange(Sex, Age_Band, Duration) %>%
  # Pivot (Wide Format: Year_1, Year_2 ...)
  tidyr::pivot_wider(
    names_from = Duration,
    values_from = Exceed_Prob,
    names_prefix = "Year_"
  )

# ------------------------------------------------------------------------------
# 4. 결과 저장
# ------------------------------------------------------------------------------
if (!dir.exists("Lapse_Table")) dir.create("Lapse_Table")

out_name <- file.path("Lapse_Table", paste0("Final_Lapse_Rates_Method_", SELECTED_METHOD, ".csv"))
write.csv(final_wide, out_name, row.names = FALSE)

cat("=======================================================\n")
cat(" [SUCCESS] 최종 탈퇴율 테이블 생성 완료\n")
cat(" 저장 경로:", out_name, "\n")
cat("=======================================================\n")
print(head(final_wide))
