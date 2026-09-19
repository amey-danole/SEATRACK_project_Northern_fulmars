# ============================================================
# Figure 5: PERS by Northern Fulmar colony
# ============================================================

# ------------------------------------------------------------
# Load packages----
# ------------------------------------------------------------

library(tidyverse)
library(ggplot2)
library(viridis)
library(scales)

# ------------------------------------------------------------
# Load data----
# ------------------------------------------------------------

input_dir <- "/Users/ameydanole/Desktop/ENS_Rennes/Rennes/argh/Amey_Danole_MS_Thesis/First_hyp/latest_right_attempt/outputs/csv/"

hyp1_df <- read.csv(
  file.path(
    input_dir,
    "correct_exposure_scores_by_individual_according_to_tracking_yr.csv"
  )
)

# Standardise colony names
hyp1_df$population <- gsub(
  "\\.",
  " ",
  hyp1_df$population
)

# Combine the two Icelandic colonies
hyp1_df$population[
  hyp1_df$population == "Combined"
] <- "Iceland"

# ------------------------------------------------------------
# Colony order----
# ------------------------------------------------------------

colony_order <- c(
  "Alkefjellet",
  "Bjørnøya",
  "Jan Mayen",
  "Iceland",
  "Faroe Islands",
  "Jarsteinen",
  "Eynhallow",
  "Isle of Canna",
  "Inishkea",
  "Little Saltee"
)

# ------------------------------------------------------------
# Include all colonies----
# ------------------------------------------------------------

plot_df <- hyp1_df %>%
  filter(
    population %in% colony_order
  ) %>%
  mutate(
    population = factor(
      population,
      levels = colony_order
    )
  )

# ------------------------------------------------------------
# Colony colour palette----
# ------------------------------------------------------------

colony_colours <- c(
  "Alkefjellet" = "#E69F00",
  "Bjørnøya" = "#F0E442",
  "Jan Mayen" = "#7CB342",
  "Iceland" = "#B8C95A",
  "Faroe Islands" = "#00A6A6",
  "Jarsteinen" = "#56B4E9",
  "Eynhallow" = "#0072B2",
  "Isle of Canna" = "#5C4FA3",
  "Inishkea" = "#8C4BB8",
  "Little Saltee" = "#CC79A7"
)

# Keep only colours used in this figure
colony_colours_h1 <- colony_colours[
  names(colony_colours) %in% colony_order
]

# ------------------------------------------------------------
# Calculate descriptive statistics----
# ------------------------------------------------------------

summary_df <- plot_df %>%
  group_by(population) %>%
  summarise(
    n = n(),
    median = median(
      exposure_score,
      na.rm = TRUE
    ),
    lower = quantile(
      exposure_score,
      0.025,
      na.rm = TRUE
    ),
    upper = quantile(
      exposure_score,
      0.975,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

write.csv(
  summary_df,
  file.path(
    input_dir,
    "Figure5_colony_PERS_summary.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# Figure 5----
# ------------------------------------------------------------

figure5 <- ggplot(
  plot_df,
  aes(
    x = population,
    y = exposure_score,
    fill = population
  )
) +
  
  # Individual observations
  geom_jitter(
    aes(
      colour = population
    ),
    width = 0.16,
    height = 0,
    size = 1.8,
    alpha = 0.55
  ) +
  
  # Boxplots
  geom_boxplot(
    width = 0.62,
    colour = "black",
    linewidth = 0.5,
    alpha = 0.85,
    outlier.shape = NA
  ) +
  
  # Colony colours
  scale_fill_manual(
    values = colony_colours_h1,
    breaks = colony_order,
    drop = FALSE
  ) +
  
  scale_colour_manual(
    values = colony_colours_h1,
    breaks = colony_order,
    drop = FALSE
  ) +
  
  # Y-axis
  scale_y_continuous(
    name = "PERS",
    expand = expansion(
      mult = c(0.02, 0.08)
    )
  ) +
  
  # X-axis
  scale_x_discrete(
    name = "Colony"
  ) +
  
  # Labels
  labs(
    x = "Colony",
    y = "PERS"
  ) +
  
  # Publication theme
  theme_minimal(
    base_size = 10
  ) +
  
  theme(
    # Axes and grid
    axis.line.x = element_line(
      colour = "black",
      linewidth = 0.6
    ),
    axis.line.y = element_line(
      colour = "black",
      linewidth = 0.6
    ),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(
      colour = "grey85",
      linewidth = 0.25
    ),
    
    # Axis titles
    axis.title.x = element_text(
      size = 16,
      face = "bold",
      colour = "black",
      margin = margin(
        t = 12
      )
    ),
    axis.title.y = element_text(
      size = 16,
      face = "bold",
      colour = "black",
      margin = margin(
        r = 12
      )
    ),
    
    # Axis text
    axis.text.x = element_text(
      size = 14,
      face = "bold",
      colour = "black",
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    axis.text.y = element_text(
      size = 14,
      face = "bold",
      colour = "black"
    ),
    
    # Legend
    legend.position = "none",
    
    # Plot margins
    plot.margin = margin(
      10,
      15,
      10,
      15
    )
  )

# Display
print(figure5)

# ------------------------------------------------------------
# Save figure----
# ------------------------------------------------------------

output_file <- "/Users/ameydanole/Desktop/ENS_Rennes/Publication/First_hypothesis/Outputs/Figures/Figure_5_PERS_by_colony.png"

ggplot2::ggsave(
  filename = output_file,
  plot = figure6,
  width = 12,
  height = 8,
  units = "in",
  dpi = 900,
  bg = "white"
)
