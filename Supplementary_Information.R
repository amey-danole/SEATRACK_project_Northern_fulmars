# ============================================================
# Supplementary Information
# Tracking effort, sample sizes and model summaries
# ============================================================

# ------------------------------------------------------------
# Packages----
# ------------------------------------------------------------

library(dplyr)
library(lubridate)
library(readr)


# ------------------------------------------------------------
# File paths----
# ------------------------------------------------------------

data_dir <- "/Users/ameydanole/Desktop/ENS_Rennes/Rennes/argh/Amey_Danole_MS_Thesis/Month/input_data"

output_dir <- "/Users/ameydanole/Desktop/ENS_Rennes/Publication/Supplementary_Information/csv"


# ------------------------------------------------------------
# Load tracking data----
# ------------------------------------------------------------

tracking_data <- readRDS(
  file.path(
    data_dir,
    "SEATRACK_FUGLA_20220307_v2.3_FA.rds"
  )
)

summary_info <- readRDS(
  file.path(
    data_dir,
    "summaryTable.rds"
  )
)


# ------------------------------------------------------------
# Merge tracking and colony information----
# ------------------------------------------------------------

input_df <- tracking_data %>%
  left_join(
    summary_info,
    by = "ring"
  ) %>%
  rename(
    individ_id = ring
  ) %>%
  mutate(
    year = year(timestamp),
    month = month(timestamp),
    tracking_year = if_else(
      month < 7,
      year - 1L,
      year
    )
  )


# ------------------------------------------------------------
# Tracking effort by colony and tracking year----
# Full tracking dataset
# ------------------------------------------------------------

tracking_year_summary <- input_df %>%
  group_by(
    colony,
    tracking_year
  ) %>%
  summarise(
    number_of_birds_tracked = n_distinct(individ_id),
    number_of_relocations = n(),
    .groups = "drop"
  ) %>%
  arrange(
    colony,
    tracking_year
  )


# ------------------------------------------------------------
# Tracking effort by colony and calendar month----
# Full tracking dataset
# ------------------------------------------------------------

month_summary <- input_df %>%
  group_by(
    colony,
    month
  ) %>%
  summarise(
    number_of_birds_tracked = n_distinct(individ_id),
    number_of_relocations = n(),
    .groups = "drop"
  ) %>%
  arrange(
    colony,
    month
  )


# ------------------------------------------------------------
# Tracking effort by individual----
# Full tracking dataset
# ------------------------------------------------------------

individual_summary <- input_df %>%
  group_by(
    colony,
    individ_id
  ) %>%
  summarise(
    number_of_relocations = n(),
    .groups = "drop"
  ) %>%
  arrange(
    colony,
    individ_id
  )


# ------------------------------------------------------------
# Restrict to the non-breeding period used in H1----
# October-March
# ------------------------------------------------------------

non_breeding_df <- input_df %>%
  filter(
    month %in% c(10, 11, 12, 1, 2, 3)
  )


# ------------------------------------------------------------
# Tracking effort by colony and tracking year----
# Non-breeding period used in H1
# ------------------------------------------------------------

tracking_year_non_breeding_summary <- non_breeding_df %>%
  group_by(
    colony,
    tracking_year
  ) %>%
  summarise(
    number_of_birds_tracked = n_distinct(individ_id),
    number_of_relocations = n(),
    .groups = "drop"
  ) %>%
  arrange(
    colony,
    tracking_year
  )


# ------------------------------------------------------------
# Create output directory if necessary----
# ------------------------------------------------------------

if (!dir.exists(output_dir)) {
  dir.create(
    output_dir,
    recursive = TRUE
  )
}


# ------------------------------------------------------------
# Write supplementary tables----
# ------------------------------------------------------------

write_csv(
  tracking_year_summary,
  file.path(
    output_dir,
    "tracking_effort_by_colony_tracking_year.csv"
  )
)

write_csv(
  month_summary,
  file.path(
    output_dir,
    "tracking_effort_by_colony_month.csv"
  )
)

write_csv(
  individual_summary,
  file.path(
    output_dir,
    "tracking_effort_by_individual.csv"
  )
)

write_csv(
  tracking_year_non_breeding_summary,
  file.path(
    output_dir,
    "tracking_effort_by_colony_tracking_year_non_breeding.csv"
  )
)


# ============================================================
# Print supplementary tables for inspection----
# ============================================================

tracking_year_summary <- read_csv(
  file.path(
    output_dir,
    "tracking_effort_by_colony_tracking_year.csv"
  ),
  show_col_types = FALSE
)

month_summary <- read_csv(
  file.path(
    output_dir,
    "tracking_effort_by_colony_month.csv"
  ),
  show_col_types = FALSE
)

individual_summary <- read_csv(
  file.path(
    output_dir,
    "tracking_effort_by_individual.csv"
  ),
  show_col_types = FALSE
)

tracking_year_non_breeding_summary <- read_csv(
  file.path(
    output_dir,
    "tracking_effort_by_colony_tracking_year_non_breeding.csv"
  ),
  show_col_types = FALSE
)


# ------------------------------------------------------------
# Tracking year table----
# Full tracking dataset
# ------------------------------------------------------------

cat("\n\n===== TRACKING YEAR SUMMARY =====\n\n")

print(
  tracking_year_summary,
  n = Inf
)


# ------------------------------------------------------------
# Month table----
# Full tracking dataset
# ------------------------------------------------------------

cat("\n\n===== MONTH SUMMARY =====\n\n")

print(
  month_summary,
  n = Inf
)


# ------------------------------------------------------------
# Individual table----
# Full tracking dataset
# ------------------------------------------------------------

cat("\n\n===== INDIVIDUAL SUMMARY =====\n\n")

print(
  individual_summary,
  n = Inf
)


# ------------------------------------------------------------
# Non-breeding tracking year table----
# October-March, used in H1
# ------------------------------------------------------------

cat("\n\n===== NON-BREEDING TRACKING YEAR SUMMARY =====\n\n")

print(
  tracking_year_non_breeding_summary,
  n = Inf
)

