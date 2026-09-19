# ============================================================
# Figure 6: PERS by PAME region
# ============================================================

# ------------------------------------------------------------
# Load packages----
# ------------------------------------------------------------

library(tidyverse)
library(ggplot2)
library(scales)

# ------------------------------------------------------------
# Load monthly PERS data----
# ------------------------------------------------------------

input_dir <- "/Users/ameydanole/Desktop/ENS_Rennes/Publication/Second_hypothesis/Outputs/Tables/csv/"

pers_df <- read.csv(
  file.path(
    input_dir,
    "PERS_by_region_and_month.csv"
  )
)

# ------------------------------------------------------------
# Region order----
# ------------------------------------------------------------

region_order <- c(
  "North Sea",
  "Faroe Plateau",
  "Iceland Shelf and Sea",
  "Canadian Eastern Arctic - West Greenland",
  "Barents Sea"
)

# ------------------------------------------------------------
# Region colour palette----
# ------------------------------------------------------------

region_colours <- c(
  "North Sea" = "#E69F00",
  "Faroe Plateau" = "#F0E442",
  "Iceland Shelf and Sea" = "#7CB342",
  "Canadian Eastern Arctic - West Greenland" = "#00A6A6",
  "Barents Sea" = "#0072B2"
)

# ------------------------------------------------------------
# Prepare plotting data----
# ------------------------------------------------------------

plot_df <- pers_df %>%
  filter(
    Region %in% region_order
  ) %>%
  mutate(
    Region = factor(
      Region,
      levels = region_order
    )
  )

# ------------------------------------------------------------
# Calculate descriptive statistics----
# ------------------------------------------------------------

summary_df <- plot_df %>%
  group_by(
    Region
  ) %>%
  summarise(
    n = n(),
    median = median(
      PERS,
      na.rm = TRUE
    ),
    lower = quantile(
      PERS,
      0.025,
      na.rm = TRUE
    ),
    upper = quantile(
      PERS,
      0.975,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

write.csv(
  summary_df,
  file.path(
    input_dir,
    "Figure6_PERS_by_PAME_region_summary.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# Figure 6----
# ------------------------------------------------------------

figure6 <- ggplot(
  plot_df,
  aes(
    x = Region,
    y = PERS,
    fill = Region
  )
) +
  
  # Individual monthly observations
  geom_jitter(
    aes(
      colour = Region
    ),
    width = 0.16,
    height = 0,
    size = 1.8,
    alpha = 0.55
  ) +
  
  # Boxplots with standard whiskers
  geom_boxplot(
    width = 0.62,
    colour = "black",
    linewidth = 0.5,
    alpha = 0.85,
    outlier.shape = NA
  ) +
  
  # Region colours
  scale_fill_manual(
    values = region_colours,
    breaks = region_order,
    drop = FALSE
  ) +
  
  scale_colour_manual(
    values = region_colours,
    breaks = region_order,
    drop = FALSE
  ) +
  
  # Region abbreviations
  scale_x_discrete(
    name = "PAME Region",
    labels = c(
      "North Sea" = "NS",
      "Faroe Plateau" = "FP",
      "Iceland Shelf and Sea" = "ISS",
      "Canadian Eastern Arctic - West Greenland" = "CEAWG",
      "Barents Sea" = "BS"
    )
  ) +
  
  # Y-axis
  scale_y_continuous(
    name = "PERS",
    expand = expansion(
      mult = c(0.02, 0.08)
    )
  ) +
  
  # Visually start the y-axis at 15 without removing observations
  coord_cartesian(
    ylim = c(15, NA)
  ) +
  
  # Labels
  labs(
    x = "PAME Region",
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
      size = 18,
      face = "bold",
      colour = "black",
      margin = margin(
        t = 12
      )
    ),
    axis.title.y = element_text(
      size = 18,
      face = "bold",
      colour = "black",
      margin = margin(
        r = 12
      )
    ),
    
    # Axis text
    axis.text.x = element_text(
      size = 16,
      face = "bold",
      colour = "black",
      angle = 0,
      hjust = 0.5,
      vjust = 0.5
    ),
    axis.text.y = element_text(
      size = 16,
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

# ------------------------------------------------------------
# Display
# ------------------------------------------------------------

print(figure6)

# ------------------------------------------------------------
# Save figure
# ------------------------------------------------------------

output_file <- "/Users/ameydanole/Desktop/ENS_Rennes/Publication/Second_hypothesis/Outputs/Figures/Figure_8_PERS_by_PAME_region.png"

ggplot2::ggsave(
  filename = output_file,
  plot = figure6,
  width = 12,
  height = 8,
  units = "in",
  dpi = 900,
  bg = "white"
)
