# Hypothesis 1 analysis

# Hypothesis 1:
# PERS varies among Northern fulmar colonies and among
# individuals within colonies during the non-breeding season.


# ============================================================
# 1. PACKAGES----
# ============================================================

library(dplyr)
library(ggplot2)
library(lme4)
library(emmeans)


# ============================================================
# 2. LOAD DATA----
# ============================================================

setwd(
  "/Users/ameydanole/Desktop/ENS_Rennes/Rennes/argh/Amey_Danole_MS_Thesis/First_hyp/latest_right_attempt/outputs/csv"
)

# Tracking-year-specific PERS values
df <- read.csv(
  "correct_exposure_scores_by_individual_according_to_tracking_yr.csv"
)

dim(df)

# Established colony-level descriptive statistics used in
# the manuscript
old_summary <- read.csv(
  "correct_ind_pers_by_population.csv"
)


# ============================================================
# 3. INSPECT INPUT DATA----
# ============================================================

str(df)
dim(df)
names(df)
head(df)

old_summary


sample_size <- read.csv(
  "/Users/ameydanole/Desktop/ENS_Rennes/Publication/First_hypothesis/Outputs/Sample_size_KDE_hyp_1.csv"
)

sample_size <- sample_size %>%
  mutate(
    population = case_when(
      colony == "Faroe Islands" ~ "Faroe.Islands",
      colony == "Jan Mayen" ~ "Jan.Mayen",
      colony == "Isle of Canna" ~ "Isle.of.Canna",
      colony == "Iceland" ~ "Combined",
      colony == "Little Saltee" ~ "Little.Saltee",
      TRUE ~ colony
    ),
    individual = gsub("-", ".", individ_id)
  )

corrected_df <- df %>%
  left_join(
    sample_size %>%
      select(
        population,
        individual,
        tracking_year,
        Number_of_relocations
      ),
    by = c(
      "population",
      "individual",
      "tracking_year"
    )
  )

cat(
  "Missing relocation counts:",
  sum(is.na(corrected_df$Number_of_relocations)),
  "\n"
)

corrected_df <- corrected_df %>%
  filter(Number_of_relocations >= 50)

analysis_df <- corrected_df[, c(
  "population",
  "individual",
  "tracking_year",
  "exposure_score"
)]

colnames(analysis_df) <- c(
  "Colony",
  "Individ_id",
  "Tracking_year",
  "pers"
)

analysis_df_2 <- analysis_df[
  analysis_df$Colony != "Alkefjellet",
]

# ============================================================
# 6. DESCRIPTIVE STATISTICS----
# ============================================================

# Number of unique individuals represented in the H1 dataset.

individual_counts <- analysis_df %>%
  group_by(Colony) %>%
  summarise(
    n_individuals = n_distinct(Individ_id),
    .groups = "drop"
  ) %>%
  arrange(desc(n_individuals))

individual_counts


# Number of tracking-year-specific PERS observations and
# descriptive statistics calculated directly from analysis_df.

# These are provided as an audit of the underlying dataset.
# The established manuscript-level descriptive values are
# retained separately in old_summary.

pers_by_colony_audit <- analysis_df %>%
  group_by(Colony) %>%
  summarise(
    n = n(),
    mean_PERS = mean(pers, na.rm = TRUE),
    median_PERS = median(pers, na.rm = TRUE),
    variance_PERS = var(pers, na.rm = TRUE),
    SD_PERS = sd(pers, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_PERS))

pers_by_colony_audit


# Established manuscript-level colony summaries.

manuscript_pers_summary <- old_summary %>%
  select(
    population,
    pers,
    variance_pers,
    median_intra_ind_variance_pers,
    percent_nas
  )

manuscript_pers_summary


# ============================================================
# 7. CHECK INDIVIDUAL × TRACKING-YEAR SAMPLE SIZES----
# ============================================================

# Identify individual × tracking-year combinations with
# fewer than 50 relocations.

# This section requires df_mod, the processed movement dataset,
# to already be present in the R session.

if (exists("df_mod")) {
  
  sample_size_df <- df_mod %>%
    group_by(colony, individ_id, tracking_year) %>%
    summarise(
      n_tracks = n(),
      .groups = "drop"
    )
  
  sample_size_df %>%
    filter(n_tracks < 50)
  
}


# ============================================================
# 8. INTER-TRACKING-YEAR VARIATION----
# ============================================================

# Calculate variation in PERS across tracking years for each
# individual.

pers_by_individual <- analysis_df %>%
  group_by(Colony, Individ_id) %>%
  summarise(
    n_tracking_years = n(),
    individual_mean_PERS = mean(
      pers,
      na.rm = TRUE
    ),
    individual_variance_PERS = var(
      pers,
      na.rm = TRUE
    ),
    .groups = "drop"
  )


# Summarise inter-tracking-year variance by colony.

# Only individuals represented in more than one tracking year
# can contribute to this measure.

inter_year_variance <- pers_by_individual %>%
  filter(n_tracking_years > 1) %>%
  group_by(Colony) %>%
  summarise(
    mean_inter_year_variance =
      mean(
        individual_variance_PERS,
        na.rm = TRUE
      ),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_inter_year_variance))

inter_year_variance


# Established inter-tracking-year variance values reported
# in the manuscript.

manuscript_inter_year_variance <- old_summary %>%
  select(
    population,
    median_intra_ind_variance_pers
  ) %>%
  arrange(desc(median_intra_ind_variance_pers))

manuscript_inter_year_variance


# ============================================================
# 9. EXPLORATORY DATA ANALYSIS----
# ============================================================

# Distribution of PERS.

dotchart(
  analysis_df$pers,
  main = "Tracking-year-specific PERS"
)


# Potential low-PERS observations.

analysis_df %>%
  filter(pers < 15)


# PERS distribution by colony.

ggplot(
  analysis_df,
  aes(
    x = Colony,
    y = pers
  )
) +
  geom_boxplot()


# Overall PERS distribution.

ggplot(
  analysis_df,
  aes(x = pers)
) +
  geom_histogram(
    binwidth = 0.2
  )


# ============================================================
# 10. EXCLUDE ALKEFJELLET FROM THE FINAL MODEL----
# ============================================================

# Alkefjellet showed substantially greater variation in PERS
# among individuals and was associated with heteroskedasticity
# and poorer model fit.

analysis_df_2 <- analysis_df[
  analysis_df$Colony != "Alkefjellet",
]

dim(analysis_df_2)

length(
  unique(
    analysis_df_2$Individ_id
  )
)


# Confirm colonies included in the final model.

unique(
  analysis_df_2$Colony
)


# Descriptive statistics for the final modelling dataset.

analysis_df_2 %>%
  group_by(Colony) %>%
  summarise(
    n = n(),
    n_individuals = n_distinct(Individ_id),
    mean_PERS = mean(
      pers,
      na.rm = TRUE
    ),
    median_PERS = median(
      pers,
      na.rm = TRUE
    ),
    variance_PERS = var(
      pers,
      na.rm = TRUE
    ),
    SD_PERS = sd(
      pers,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_PERS))


# ============================================================
# 11. MODEL SELECTION: DISTRIBUTION AND LINK FUNCTION----
# ============================================================

control_params <- glmerControl(
  optimizer = "bobyqa",
  optCtrl = list(
    maxfun = 100000
  )
)


# Gamma distribution with reciprocal link

g_gamma_reciprocal <- glmer(
  pers ~ Colony + (1 | Individ_id),
  data = analysis_df_2,
  family = Gamma(
    link = "inverse"
  ),
  control = control_params
)


# Gamma distribution with log link

g_gamma_log <- glmer(
  pers ~ Colony + (1 | Individ_id),
  data = analysis_df_2,
  family = Gamma(
    link = "log"
  ),
  control = control_params
)


# Gamma distribution with square-root link

g_gamma_sqrt <- glmer(
  pers ~ Colony + (1 | Individ_id),
  data = analysis_df_2,
  family = Gamma(
    link = "sqrt"
  ),
  control = control_params
)


# Gamma distribution with identity link

g_gamma_identity <- glmer(
  pers ~ Colony + (1 | Individ_id),
  data = analysis_df_2,
  family = Gamma(
    link = "identity"
  ),
  control = control_params
)


# Inverse Gaussian distribution with identity link

g_inverse_gaussian <- glmer(
  pers ~ Colony + (1 | Individ_id),
  data = analysis_df_2,
  family = inverse.gaussian(
    link = "identity"
  ),
  control = control_params
)


# Compare candidate models using AIC.

model_AIC <- AIC(
  g_gamma_reciprocal,
  g_gamma_log,
  g_gamma_sqrt,
  g_gamma_identity,
  g_inverse_gaussian
)

model_AIC


# Compare candidate models using BIC.

model_BIC <- BIC(
  g_gamma_reciprocal,
  g_gamma_log,
  g_gamma_sqrt,
  g_gamma_identity,
  g_inverse_gaussian
)

model_BIC


# ============================================================
# 12. FINAL H1 MODEL----
# ============================================================

g_mod_2 <- glmer(
  pers ~ Colony + (1 | Individ_id),
  data = analysis_df_2,
  family = Gamma(
    link = "identity"
  ),
  control = control_params
)


# Model summary.

summary(g_mod_2)


# Model information criteria.

AIC(g_mod_2)

BIC(g_mod_2)


# Random-effect variance.

VarCorr(g_mod_2)


# ============================================================
# 13. ESTIMATED MARGINAL MEANS----
# ============================================================

emm <- emmeans(
  g_mod_2,
  ~ Colony
)

summary(emm)


# ============================================================
# 14. PAIRWISE COLONY COMPARISONS----
# ============================================================

pairwise_colony_comparisons <- pairs(
  emm,
  adjust = "tukey"
)

pairwise_colony_comparisons

summary(
  pairwise_colony_comparisons
)


# ============================================================
# 15. MODEL DIAGNOSTICS----
# ============================================================

# Standard model diagnostic plots.

plot(g_mod_2)


# Residual distribution.

qqnorm(
  resid(g_mod_2)
)

qqline(
  resid(g_mod_2)
)


hist(
  resid(g_mod_2),
  breaks = 40,
  main = "Residual distribution",
  xlab = "Residuals"
)


# Residuals versus fitted values.

plot(
  fitted(g_mod_2),
  resid(g_mod_2),
  xlab = "Fitted values",
  ylab = "Residuals"
)

abline(
  h = 0
)


# Residuals by colony.

residual_df <- data.frame(
  Colony = analysis_df_2$Colony,
  residuals = resid(g_mod_2)
)

ggplot(
  residual_df,
  aes(
    x = Colony,
    y = residuals
  )
) +
  geom_boxplot()


# ============================================================
# 16. SAVE FINAL MODEL OUTPUTS----
# ============================================================

tables_dir <- "/Users/ameydanole/Desktop/ENS_Rennes/Publication/First_hypothesis/Outputs/Tables"


write.csv(
  as.data.frame(
    summary(emm)
  ),
  file.path(
    tables_dir,
    "H1_estimated_marginal_means.csv"
  ),
  row.names = FALSE
)


write.csv(
  as.data.frame(
    summary(pairwise_colony_comparisons)
  ),
  file.path(
    tables_dir,
    "H1_pairwise_colony_comparisons.csv"
  ),
  row.names = FALSE
)
