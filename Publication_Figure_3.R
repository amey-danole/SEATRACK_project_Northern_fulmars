
# ============================================================
# FIGURE 3 - GLOBAL DISTRIBUTION OF MARINE FLOATING PLASTIC
# ============================================================


# Loading essential packages----
library(spData)
library(tidyverse)
library(ggplot2)
library(leaflet)
library(rnaturalearth)
library(sf)
library(sp)
library(raster)
library(mapproj)
library(rayshader)
library(maps)
library(RColorBrewer)
library(viridis)
library(cowplot)
library(adehabitatHR)
library(dplyr)
library(terra)


# Loading data----

datadir <- "/Users/ameydanole/Desktop/ENS_Rennes/Rennes/argh/Amey_Danole_MS_Thesis/Ind/input_data"

mylocs <- readRDS(
  paste0(
    datadir,
    "/SEATRACK_FUGLA_20220307_v2.3_FA.rds"
  )
)

summary_info <- readRDS(
  paste0(
    datadir,
    "/summaryTable.rds"
  )
)

indiv_merged_df <- merge(
  mylocs,
  summary_info,
  by = "ring"
)

nbs_mylocs <- indiv_merged_df %>%
  dplyr::mutate(
    year = year(timestamp)
  )

names(nbs_mylocs)[
  names(nbs_mylocs) == "ring"
] <- "individ_id"

all_data <- nbs_mylocs

df <- all_data[
  !is.na(all_data$timestamp),
]

df <- df %>%
  dplyr::mutate(
    month = month(timestamp)
  )

months <- sort(
  unique(df$month)
)

df_mod <- df %>%
  mutate(
    tracking_year = ifelse(
      month < 7,
      year - 1,
      year
    )
  )


# Convert tracking data to sf----

df_mod_sf <- df_mod %>%
  sf::st_as_sf(
    coords = c(
      "col_lon",
      "col_lat"
    ),
    crs = 4326
  )

ind_merge_sf <- df_mod_sf


# Load world map data----

world <- rnaturalearth::ne_countries(
  scale = "medium",
  returnclass = "sf"
)

world <- sf::st_transform(
  world,
  "+proj=longlat +datum=WGS84"
)


# Define land----

land <- as(
  world,
  "Spatial"
)


# Define WGS84 projection----

proj_wgs84 <- sp::CRS(
  sp::proj4string(land)
)


# Loading plastic raster----

setwd(
  datadir
)

plastics <- raster::raster(
  "00_PlasticsRaster.tif"
)


# Normalise plastic distribution----

plastics2 <- plastics

plastics2[
  is.na(plastics2)
] <- 0

p_sum1 <- plastics2 /
  sum(
    raster::getValues(
      plastics2
    )
  )

p_sum1[
  is.na(plastics)
] <- NA


# Define Robinson projection----

proj <- "+proj=robin"


# Project land to Robinson projection----

land_sf <- sf::st_as_sf(
  land
)

world_prj <- land_sf %>%
  sf::st_transform(
    proj
  )


# Define global map extent----

CP <- sf::st_bbox(
  c(
    xmin = -180,
    xmax = 180,
    ymin = -90,
    ymax = 90
  )
) %>%
  sf::st_as_sfc() %>%
  sf::st_segmentize(
    1
  ) %>%
  sf::st_set_crs(
    4326
  )

CP_prj <- CP %>%
  sf::st_transform(
    crs = proj
  )


# Get transformed bounding box and expand----

xlim <- sf::st_bbox(
  CP_prj
)[
  c(
    "xmin",
    "xmax"
  )
] * 1.2

ylim <- sf::st_bbox(
  CP_prj
)[
  c(
    "ymin",
    "ymax"
  )
] * 1.2


# Create enclosing rectangle----

encl_rect <- list(
  cbind(
    c(
      xlim[1],
      xlim[2],
      xlim[2],
      xlim[1],
      xlim[1]
    ),
    c(
      ylim[1],
      ylim[1],
      ylim[2],
      ylim[2],
      ylim[1]
    )
  )
) %>%
  sf::st_polygon() %>%
  sf::st_sfc(
    crs = proj
  )


# Calculate area outside Earth outline----

cookie <- sf::st_difference(
  encl_rect,
  CP_prj
)


# Project plastic raster to Robinson projection ----

# Convert raster to terra SpatRaster
plastics_terra <- terra::rast(plastics)

# Remove the exact polar rows
plastics_terra <- terra::crop(
  plastics_terra,
  terra::ext(
    terra::xmin(plastics_terra),
    terra::xmax(plastics_terra),
    -89.5,
    89.5
  )
)

# Project to Robinson
p_dens_proj <- terra::project(
  plastics_terra,
  proj,
  method = "bilinear"
)

# Convert raster to dataframe for ggplot
p_df <- terra::as.data.frame(
  p_dens_proj,
  xy = TRUE,
  na.rm = TRUE
) %>%
  tibble::as_tibble()

# Rename plastic-value column
names(p_df)[3] <- "plastic"

# Figure 3: Global distribution of marine floating plastic debris----

sp3 <- ggplot() +
  
  # Floating plastic debris
  geom_raster(
    data = p_df,
    aes(
      x = x,
      y = y,
      fill = plastic
    )
  ) +
  
  # Land
  geom_sf(
    data = world_prj,
    fill = "grey75",
    colour = "grey45",
    linewidth = 0.25
  ) +
  
  # Mask outside the global map
  geom_sf(
    data = cookie,
    fill = "white",
    colour = NA
  ) +
  
  # Plastic colour scale
  scale_fill_viridis_c(
    option = "inferno",
    direction = -1,
    name = "log (plastic count / km²)",
    guide = guide_colourbar(
      title.position = "top",
      title.hjust = 0.5,
      direction = "horizontal",
      barwidth = unit(8, "cm"),
      barheight = unit(0.45, "cm")
    )
  ) +
  
  # Robinson projection
  coord_sf(
    expand = FALSE
  ) +
  
  labs(
    x = NULL,
    y = NULL
  ) +
  
  theme_minimal(
    base_size = 10
  ) +
  
  theme(
    # Remove longitude/latitude axes and ticks
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    
    # Remove grid
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    
    # Remove panel border
    panel.border = element_blank(),
    
    # Legend
    legend.position = "top",
    legend.direction = "horizontal",
    legend.title = element_text(
      size = 12,
      face = "bold",
      colour = "black"
    ),
    legend.text = element_text(
      size = 10,
      face = "bold",
      colour = "black"
    ),
    legend.box.margin = margin(
      5,
      0,
      0,
      0
    ),
    
    # Plot margins
    plot.margin = margin(
      10,
      15,
      10,
      15
    )
  )

# Display figure----

print(sp3)

# ============================================================
# Save figure----
# ============================================================

output_file <- "/Users/ameydanole/Desktop/ENS_Rennes/Publication/Figures/Inferno_Global_distribution_of_marine_floating_plastic_debris.png"

ggplot2::ggsave(
  filename = output_file,
  plot = sp3,
  width = 14,
  height = 8,
  units = "in",
  dpi = 900,
  bg = "white"
)
