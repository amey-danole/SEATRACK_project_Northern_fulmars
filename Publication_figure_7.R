# ============================================================
# Figure 7: Relationship between PERS and EcoQO
# ============================================================

# ------------------------------------------------------------
# Load packages----
# ------------------------------------------------------------

library(tidyverse)
library(ggplot2)

# ------------------------------------------------------------
# Load data----
# ------------------------------------------------------------

input_file <- "/Users/ameydanole/Desktop/ENS_Rennes/Rennes/argh/Amey_Danole_MS_Thesis/Second_hyp/outputs/csv/correct_hyp_2_analysis_df.csv"

hyp2_analysis_df <- read.csv(
  input_file
)

# Standardise region names
hyp2_analysis_df$Region <- gsub(
  "\\.",
  " ",
  hyp2_analysis_df$Region
)

# ------------------------------------------------------------
# PAME region order----
# ------------------------------------------------------------

region_order <- c(
  "North Sea",
  "Faroe Plateau",
  "Iceland Shelf and Sea",
  "Canadian Eastern Arctic - West Greenland",
  "Barents Sea"
)

hyp2_analysis_df <- hyp2_analysis_df %>%
  mutate(
    Region = factor(
      Region,
      levels = region_order
    )
  )

# ------------------------------------------------------------
# PAME region colours----
# ------------------------------------------------------------

region_colours <- c(
  "North Sea" = "#E69F00",
  "Faroe Plateau" = "#F0E442",
  "Iceland Shelf and Sea" = "#7CB342",
  "Canadian Eastern Arctic - West Greenland" = "#00A6A6",
  "Barents Sea" = "#0072B2"
)

# ------------------------------------------------------------
# Linear regression----
# ------------------------------------------------------------

hyp2_lm <- lm(
  EcoQO.value ~ median_pers,
  data = hyp2_analysis_df
)

# ------------------------------------------------------------
# Figure 7----
# ------------------------------------------------------------

figure7 <- ggplot(
  hyp2_analysis_df,
  aes(
    x = median_pers,
    y = EcoQO.value,
    fill = Region
  )
) +
  
  # Linear regression line across all PAME regions
  geom_smooth(
    aes(
      x = median_pers,
      y = EcoQO.value,
      group = 1
    ),
    method = "lm",
    formula = y ~ x,
    se = FALSE,
    colour = "black",
    linewidth = 1.0,
    inherit.aes = FALSE
  ) + 
  
  # PAME region observations
  geom_point(
    aes(
      fill = Region
    ),
    shape = 21,
    size = 6,
    colour = "black",
    stroke = 0.6
  ) +
  
  # Region colours
  scale_fill_manual(
    values = region_colours,
    breaks = region_order,
    drop = FALSE
  ) +
  
  # Axis labels
  labs(
    x = "PERS",
    y = "EcoQO",
    fill = "PAME Region"
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
    panel.grid.major.x = element_line(
      colour = "grey85",
      linewidth = 0.25
    ),
    panel.grid.major.y = element_line(
      colour = "grey85",
      linewidth = 0.25
    ),
    panel.grid.minor = element_blank(),
    
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
      colour = "black"
    ),
    axis.text.y = element_text(
      size = 14,
      face = "bold",
      colour = "black"
    ),
    
    # Legend
    legend.position = "right",
    legend.title = element_text(
      size = 16,
      face = "bold",
      colour = "black"
    ),
    legend.text = element_text(
      size = 14,
      face = "bold",
      colour = "black"
    ),
    legend.key.height = unit(
      0.8,
      "cm"
    ),
    legend.key.width = unit(
      0.8,
      "cm"
    ),
    
    # Plot margins
    plot.margin = margin(
      10,
      15,
      10,
      10
    )
  )

# ------------------------------------------------------------
# Display----
# ------------------------------------------------------------

print(figure7)

# ------------------------------------------------------------
# Save figure----
# ------------------------------------------------------------

output_file <- "/Users/ameydanole/Desktop/ENS_Rennes/Publication/First_hypothesis/Outputs/Figures/Figure_7_EcoQO_vs_PERS.png"

ggplot2::ggsave(
  filename = output_file,
  plot = figure7,
  width = 12,
  height = 8,
  units = "in",
  dpi = 900,
  bg = "white"
)
