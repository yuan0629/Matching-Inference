# Run this file from its own folder. It creates the single analysis-ready CSV.
file_arg <- grep("^--file=", commandArgs(), value = TRUE)
if (length(file_arg)) {
  setwd(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
}

raw <- read.csv(
  "raw data.csv", stringsAsFactors = FALSE,
  na.strings = c("", "NA"), check.names = FALSE
)
raw$date <- as.Date(raw$date)
as_bool <- function(x) tolower(as.character(x)) == "true"

# Historical outcomes used only for the 2015--2024 pair-level comparison.
historical <- raw[
  raw$season %in% 2015:2024 & raw$game_type == "R" &
    !is.na(raw$home_team) & !is.na(raw$umpire) &
    !is.na(raw$accuracy_above_expected),
]
historical$sample <- "historical_2015_2024"
historical$time_index <- NA_integer_

# Main 2025 panel: regular season, fixed Opening Day umpire pool, observed AAE,
# and no daily row or column duplication. Doubleheader conflicts were removed
# by the primary_2025_keep flag when the raw source table was assembled.
main <- raw[
  raw$season == 2025 & raw$game_type == "R" &
    as_bool(raw$primary_2025_keep) &
    !is.na(raw$accuracy_above_expected) &
    !is.na(raw$row_index) & !is.na(raw$col_index),
]
main <- main[order(main$date, main$row_index), ]
main$sample <- "analysis_2025"
main$time_index <- match(main$date, sort(unique(main$date)))

stopifnot(
  nrow(historical) == 22639L,
  nrow(main) == 2073L,
  length(unique(main$time_index)) == 180L,
  length(unique(main$row_index)) == 30L,
  length(unique(main$col_index)) == 75L,
  !anyDuplicated(main[c("date", "row_index")]),
  !anyDuplicated(main[c("date", "col_index")])
)

keep <- c(
  "sample", "date", "time_index", "game_pk", "game_number",
  "is_doubleheader", "row_index", "col_index", "home_team",
  "home_team_name", "away_team", "away_team_name", "umpire",
  "called_pitches", "correct_calls", "incorrect_calls",
  "correct_calls_above_expected", "accuracy", "expected_accuracy",
  "accuracy_above_expected", "consistency", "favor_home",
  "total_run_impact", "home_batter_impact", "home_pitcher_impact",
  "away_batter_impact", "away_pitcher_impact", "source_url"
)
clean <- rbind(historical[keep], main[keep])
write.csv(clean, "clean data.csv", row.names = FALSE, na = "")

cat("Historical observations:", nrow(historical), "\n")
cat("2025 analysis observations:", nrow(main), "\n")
cat("2025 dimensions: d1 = 30, d2 = 75, T = 180\n")
