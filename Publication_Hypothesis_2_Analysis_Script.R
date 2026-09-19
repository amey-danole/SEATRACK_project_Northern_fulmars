# ============================================================
# HYPOTHESIS 2
# ============================================================
#
# Hypothesis 2:
# Movement-derived Plastic Exposure Risk Scores (PERS) are
# associated with established EcoQO benchmarks across marine regions.
#
# Analytical workflow:
# 1. Assign SEATRACK relocations to PAME/LME regions.
# 2. Calculate monthly Utilisation Distributions (UDs) for
#    each region using kernel density estimation.
# 3. Extract the 95% UDs.
# 4. Overlay each regional UD with the global floating-plastic
#    raster to calculate a monthly PERS.
# 5. Quantify the proportion of the regional UD for which
#    plastic data are unavailable (%NA).
# 6. Calculate median monthly PERS and %NA for each region.
# 7. Compare regional PERS with established EcoQO values
#    using Kendall's rank correlation.
#
# ============================================================


# ============================================================
# 1. PACKAGES----
# ============================================================

library(dplyr)
library(sf)
library(lubridate)
library(rnaturalearth)
library(raster)
library(sp)
library(adehabitatHR)
library(geosphere)


# ============================================================
# 2. DIRECTORIES----
# ============================================================

data_dir <- "/Users/ameydanole/Desktop/ENS_Rennes/Rennes/argh/Amey_Danole_MS_Thesis"

output_dir <- "/Users/ameydanole/Desktop/ENS_Rennes/Publication/Second_hypothesis/Outputs/Tables"

output_csv_dir <- file.path(
  output_dir,
  "csv"
)

output_kde_dir <- file.path(
  output_dir,
  "unique_tifs"
)

output_map_dir <- file.path(
  output_dir,
  "unique_distributions"
)

output_na_map_dir <- file.path(
  output_dir,
  "na_maps"
)

plastic_dir <- file.path(
  data_dir,
  "Month",
  "input_data"
)

dir.create(output_csv_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_kde_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_map_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_na_map_dir, recursive = TRUE, showWarnings = FALSE)


# ============================================================
# 3. LOAD SEATRACK DATA----
# ============================================================

mylocs <- readRDS(
  file.path(
    data_dir,
    "Ind",
    "input_data",
    "SEATRACK_FUGLA_20220307_v2.3_FA.rds"
  )
)

summary_info <- readRDS(
  file.path(
    data_dir,
    "Ind",
    "input_data",
    "summaryTable.rds"
  )
)


# ============================================================
# 4. MERGE AND CURATE TRACKING DATA----
# ============================================================

indiv_merged_df <- merge(
  mylocs,
  summary_info,
  by = "ring"
) %>%
  mutate(
    year = year(timestamp)
  )

names(indiv_merged_df)[
  names(indiv_merged_df) == "ring"
] <- "individ_id"


# Combine the Skjalfandi and Langanes colonies as Iceland.
# The two colonies are approximately 130 km apart and were
# therefore treated as a single colony in the analysis.

colonies_to_combine <- c(
  "Skjalfandi",
  "Langanes"
)

indiv_merged_df$colony[
  indiv_merged_df$colony %in% colonies_to_combine
] <- "Iceland"


# ============================================================
# 5. LOAD PAME REGIONS----
# ============================================================

PAME_shapefile <- st_read(
  file.path(
    data_dir,
    "Ind",
    "input_data",
    "PAME",
    "modified_LME.shp"
  ),
  quiet = TRUE
)

# Repair invalid polygon geometries.

valid_pame <- st_make_valid(
  PAME_shapefile
)


# ============================================================
# 6. ASSIGN TRACKING RELOCATIONS TO PAME REGIONS----
# ============================================================

# Convert tracking locations to sf points.

tracking_sf <- indiv_merged_df %>%
  st_as_sf(
    coords = c(
      "lon",
      "lat"
    ),
    crs = 4326
  )


# Spatially join relocations to PAME/LME polygons.

joined <- st_join(
  tracking_sf,
  valid_pame,
  join = st_intersects
)


# Retain only variables required for the regional analysis.

region_assignments <- joined %>%
  st_drop_geometry() %>%
  select(
    individ_id,
    timestamp,
    LME_NAME
  )


# Save regional assignments.

write.csv(
  region_assignments,
  file.path(
    output_csv_dir,
    "relocations_with_pame_regions.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 7. DEFINE REGIONS USED FOR HYPOTHESIS 2----
# ============================================================

relevant_regions <- c(
  "North Sea",
  "Faroe Plateau",
  "Iceland Shelf and Sea",
  "Canadian Eastern Arctic - West Greenland",
  "Barents Sea"
)


# ============================================================
# 8. PREPARE MONTHLY TRACKING DATA----
# ============================================================

# Reload the tracking data used for the monthly analysis.

monthly_mylocs <- readRDS(
  file.path(
    data_dir,
    "Month",
    "input_data",
    "SEATRACK_FUGLA_20220307_v2.3_FA.rds"
  )
)

monthly_summary <- readRDS(
  file.path(
    data_dir,
    "Month",
    "input_data",
    "summaryTable.rds"
  )
)


monthly_data <- monthly_mylocs %>%
  left_join(
    monthly_summary,
    by = "ring"
  ) %>%
  mutate(
    year = year(timestamp),
    month = month(timestamp)
  )

names(monthly_data)[
  names(monthly_data) == "ring"
] <- "individ_id"


# Add PAME region assignments.

monthly_data <- monthly_data %>%
  left_join(
    region_assignments,
    by = c(
      "individ_id",
      "timestamp"
    )
  )


# Retain only relocations assigned to regions relevant to H2.

monthly_data <- monthly_data %>%
  filter(
    LME_NAME %in% relevant_regions
  )


# ============================================================
# 9. DEFINE KERNEL GRID AND PROJECTION----
# ============================================================

land <- ne_countries(
  scale = "medium",
  returnclass = "sf"
)

land <- as(
  land,
  "Spatial"
)

proj_wgs84 <- CRS(
  proj4string(land)
)


# ============================================================
# 10. CALCULATE MONTHLY REGIONAL UTILISATION DISTRIBUTIONS----
# ============================================================

kernel_sample_sizes <- data.frame(
  Region = character(),
  Month = integer(),
  Sample_size_for_kde = integer()
)


months <- sort(
  unique(
    monthly_data$month
  )
)


for (region in relevant_regions) {
  
  region_data <- monthly_data %>%
    filter(
      LME_NAME == region
    )
  
  for (month_number in months) {
    
    tracks_month <- region_data %>%
      filter(
        month == month_number
      )
    
    n_tracks <- nrow(
      tracks_month
    )
    
    
    # Record sample size.
    
    kernel_sample_sizes <- rbind(
      kernel_sample_sizes,
      data.frame(
        Region = region,
        Month = month_number,
        Sample_size_for_kde = n_tracks
      )
    )
    
    
    # Calculate a KDE only when more than 50 relocations
    # are available for the region-month combination.
    
    if (n_tracks <= 50) {
      next
    }
    
    
    # --------------------------------------------------------
    # Define spatial extent
    # --------------------------------------------------------
    
    lon_min <- ifelse(
      min(tracks_month$lon) <= -179,
      -180,
      floor(min(tracks_month$lon)) - 1
    )
    
    lon_max <- ifelse(
      max(tracks_month$lon) >= 179,
      180,
      ceiling(max(tracks_month$lon)) + 1
    )
    
    lat_min <- ifelse(
      min(tracks_month$lat) <= -89,
      -90,
      floor(min(tracks_month$lat)) - 1
    )
    
    lat_max <- ifelse(
      max(tracks_month$lat) >= 89,
      90,
      ceiling(max(tracks_month$lat)) + 1
    )
    
    
    # --------------------------------------------------------
    # Create initial grid in WGS84
    # --------------------------------------------------------
    
    so_grid <- expand.grid(
      LON = seq(
        lon_min,
        lon_max,
        by = 1
      ),
      LAT = seq(
        lat_min,
        lat_max,
        by = 1
      )
    )
    
    coordinates(
      so_grid
    ) <- ~ LON + LAT
    
    crs(
      so_grid
    ) <- proj_wgs84
    
    
    # --------------------------------------------------------
    # Create Lambert azimuthal equal-area projection
    # centred on the regional monthly distribution
    # --------------------------------------------------------
    
    mean_loc <- geosphere::geomean(
      cbind(
        tracks_month$lon,
        tracks_month$lat
      )
    )
    
    DgProj <- CRS(
      paste0(
        "+proj=laea +lon_0=",
        mean_loc[1],
        " +lat_0=",
        mean_loc[2]
      )
    )
    
    
    # --------------------------------------------------------
    # Create 10-km kernel grid
    # --------------------------------------------------------
    
    so_grid_proj <- spTransform(
      so_grid,
      CRS = DgProj
    )
    
    coords <- so_grid_proj@coords
    
    x_min <- min(coords[, 1]) - 1000000
    x_max <- max(coords[, 1]) + 1000000
    
    y_min <- min(coords[, 2]) - 1000000
    y_max <- max(coords[, 2]) + 1000000
    
    x_seq <- seq(
      x_min,
      x_max,
      by = 10000
    )
    
    y_seq <- seq(
      y_min,
      y_max,
      by = 10000
    )
    
    null_grid <- expand.grid(
      x = x_seq,
      y = y_seq
    )
    
    coordinates(
      null_grid
    ) <- ~ x + y
    
    gridded(
      null_grid
    ) <- TRUE
    
    
    # --------------------------------------------------------
    # Convert tracking locations to projected coordinates
    # --------------------------------------------------------
    
    coordinates(
      tracks_month
    ) <- ~ lon + lat
    
    crs(
      tracks_month
    ) <- proj_wgs84
    
    tracks_projected <- spTransform(
      tracks_month,
      CRS = DgProj
    )
    
    tracks_projected$month <- factor(
      tracks_projected@data$month
    )
    
    
    # --------------------------------------------------------
    # Kernel density estimation
    # --------------------------------------------------------
    
    kud <- adehabitatHR::kernelUD(
      tracks_projected[, "month"],
      grid = null_grid,
      h = 200000
    )
    
    
    # --------------------------------------------------------
    # Extract 95% utilisation distribution
    # --------------------------------------------------------
    
    volume_ud <- adehabitatHR::getvolumeUD(
      kud
    )
    
    volume_raster <- volume_ud[[1]]
    
    hr95 <- as.data.frame(
      volume_raster
    )[, 1]
    
    hr95 <- as.numeric(
      hr95 <= 95
    )
    
    hr95 <- data.frame(
      hr95
    )
    
    coordinates(
      hr95
    ) <- coordinates(
      volume_raster
    )
    
    gridded(
      hr95
    ) <- TRUE
    
    
    # --------------------------------------------------------
    # Convert kernel to raster
    # --------------------------------------------------------
    
    kde_spixdf <- adehabitatHR::estUDm2spixdf(
      kud
    )
    
    kernel_raster <- raster(
      kde_spixdf
    )
    
    hr95_raster <- raster(
      hr95
    )
    
    
    # Retain the 95% utilisation distribution.
    
    regional_ud <- kernel_raster *
      hr95_raster
    
    
    # Standardise the regional utilisation distribution.
    
    regional_ud <- regional_ud /
      sum(
        getValues(
          regional_ud
        )
      )
    
    regional_ud[
      regional_ud == 0
    ] <- NA
    
    
    # --------------------------------------------------------
    # Crop raster to the non-NA region
    # --------------------------------------------------------
    
    na_matrix <- is.na(
      as.matrix(
        regional_ud
      )
    )
    
    columns_non_na <- which(
      colSums(na_matrix) != nrow(regional_ud)
    )
    
    rows_non_na <- which(
      rowSums(na_matrix) != ncol(regional_ud)
    )
    
    
    if (
      length(columns_non_na) == 0 ||
      length(rows_non_na) == 0
    ) {
      next
    }
    
    
    cropped_extent <- raster::extent(
      regional_ud,
      r1 = rows_non_na[1] - 2,
      r2 = rows_non_na[length(rows_non_na)] + 2,
      c1 = columns_non_na[1] - 2,
      c2 = columns_non_na[length(columns_non_na)] + 2
    )
    
    
    cropped <- raster::crop(
      regional_ud,
      cropped_extent
    )
    
    cropped[
      is.na(cropped)
    ] <- 0
    
    
    # --------------------------------------------------------
    # Mask terrestrial areas
    # --------------------------------------------------------
    
    land_projected <- spTransform(
      land,
      DgProj
    )
    
    land_polygons <- as(
      land_projected,
      "SpatialPolygons"
    )
    
    
    regional_ud_masked <- raster::mask(
      cropped,
      land_polygons,
      inverse = TRUE
    )
    
    regional_ud_masked[
      is.na(regional_ud_masked)
    ] <- 0
    
    
    regional_ud_masked <-
      regional_ud_masked /
      sum(
        getValues(
          regional_ud_masked
        )
      )
    
    
    regional_ud_final <- raster::mask(
      regional_ud_masked,
      land_polygons,
      inverse = TRUE
    )
    
    
    # --------------------------------------------------------
    # Reproject regional UD to WGS84
    # --------------------------------------------------------
    
    regional_ud_wgs84 <- projectRaster(
      regional_ud_final,
      crs = proj_wgs84,
      over = FALSE
    )
    
    
    KDE_ref <- paste0(
      region,
      "_",
      month_number
    )
    
    
    # Save regional monthly UD.
    
    writeRaster(
      regional_ud_wgs84,
      filename = file.path(
        output_kde_dir,
        paste0(
          KDE_ref,
          ".tif"
        )
      ),
      format = "GTiff",
      overwrite = TRUE
    )
    
    
    # Save regional distribution map.
    
    png(
      filename = file.path(
        output_map_dir,
        paste0(
          KDE_ref,
          ".png"
        )
      ),
      width = 1399,
      height = 455
    )
    
    plot(
      regional_ud_wgs84,
      main = region
    )
    
    plot(
      land,
      add = TRUE
    )
    
    dev.off()
  }
}


# Save KDE sample sizes.

write.csv(
  kernel_sample_sizes,
  file.path(
    output_csv_dir,
    "Sample_size_for_KDE.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 11. LOAD GLOBAL FLOATING-PLASTIC RASTER----
# ============================================================

plastic_raster <- raster(
  file.path(
    plastic_dir,
    "00_PlasticsRaster.tif"
  )
)


# Replace unavailable values with zero for normalisation.

plastic_raster_zero <- plastic_raster

plastic_raster_zero[
  is.na(plastic_raster_zero)
] <- 0


# Standardise plastic values.

plastic_distribution <- plastic_raster_zero /
  sum(
    getValues(
      plastic_raster_zero
    )
  )


# Restore unavailable cells.

plastic_distribution[
  is.na(plastic_raster)
] <- NA


# ============================================================
# 12. CALCULATE AREA OF PLASTIC-RASTER CELLS----
# ============================================================

plastic_resolution <- res(
  plastic_raster
)

earth_radius <- 6371007.2

latitude <- yFromRow(
  plastic_raster,
  1:nrow(plastic_raster)
)


cell_area <- (
  sin(
    pi / 180 *
      (
        latitude +
          plastic_resolution[2] / 2
      )
  ) -
    sin(
      pi / 180 *
        (
          latitude -
            plastic_resolution[2] / 2
        )
    )
) *
  (
    plastic_resolution[1] *
      pi / 180
  ) *
  earth_radius^2


area_raster <- setValues(
  plastic_raster,
  rep(
    cell_area,
    each = ncol(plastic_raster)
  )
)


# ============================================================
# 13. MULTIPLY REGIONAL UDs WITH PLASTIC DISTRIBUTION----
# ============================================================

kde_files <- list.files(
  output_kde_dir,
  full.names = TRUE,
  pattern = "\\.tif$"
)

monthly_results <- data.frame()

na_results <- data.frame()

for (file in kde_files) {
  
  cat("Processing:", basename(file), "\n")
  
  regional_ud <- raster(file)
  
  region_month <- tools::file_path_sans_ext(
    basename(file)
  )
  
  # ----------------------------------------------------------
  # Reproject regional UD to plastic-raster grid
  # ----------------------------------------------------------
  
  regional_ud_plastic <- projectRaster(
    regional_ud,
    plastic_raster,
    method = "bilinear"
  )
  
  # ----------------------------------------------------------
  # Account for cell area
  # ----------------------------------------------------------
  
  regional_ud_area <- regional_ud_plastic *
    area_raster /
    100000000
  
  regional_ud_area[
    is.na(regional_ud_area)
  ] <- 0
  
  # ----------------------------------------------------------
  # Calculate PERS
  # ----------------------------------------------------------
  
  regional_ud_area[
    is.na(plastic_raster)
  ] <- NA
  
  exposure_surface <-
    regional_ud_area *
    plastic_distribution
  
  exposure_surface[
    is.na(exposure_surface)
  ] <- 0
  
  exposure_score <-
    round(
      sum(
        getValues(exposure_surface)
      ) *
        1000000,
      4
    )
  
  # ----------------------------------------------------------
  # Calculate percentage of regional UD cells with unavailable
  # plastic data
  # ----------------------------------------------------------
  
  ud_cells <- regional_ud_plastic *
    area_raster /
    100000000
  
  ud_cells[
    is.na(ud_cells)
  ] <- 0
  
  ud_cells[
    ud_cells > 0
  ] <- 1
  
  plastic_unavailable <- calc(
    is.na(plastic_distribution),
    fun = function(x) {
      as.numeric(x)
    }
  )
  
  unavailable_cells <-
    ud_cells *
    plastic_unavailable
  
  n_ud_cells <- sum(
    getValues(ud_cells)
  )
  
  n_unavailable_cells <- sum(
    getValues(unavailable_cells)
  )
  
  percent_na <-
    n_unavailable_cells /
    n_ud_cells *
    100
  
  # ----------------------------------------------------------
  # Extract region and month
  # ----------------------------------------------------------
  
  name_split <- strsplit(
    region_month,
    "_",
    fixed = TRUE
  )[[1]]
  
  month_number <- as.integer(
    name_split[length(name_split)]
  )
  
  region <- paste(
    name_split[-length(name_split)],
    collapse = "_"
  )
  
  # ----------------------------------------------------------
  # Store results
  # ----------------------------------------------------------
  
  monthly_results <- rbind(
    monthly_results,
    data.frame(
      Region = region,
      Month = month_number,
      PERS = exposure_score
    )
  )
  
  na_results <- rbind(
    na_results,
    data.frame(
      Region = region,
      Month = month_number,
      Percent_NA = percent_na
    )
  )
}
  
  # ----------------------------------------------------------
  # Save exposure surface
  # ----------------------------------------------------------
  
  writeRaster(
    regional_ud_area,
    filename = file.path(
      output_dir,
      "multiplication_rasters",
      paste0(
        region_month,
        ".tif"
      )
    ),
    format = "GTiff",
    overwrite = TRUE
  )


# ============================================================
# 14. SAVE MONTHLY RESULTS----
# ============================================================

write.csv(
  monthly_results,
  file.path(
    output_csv_dir,
    "PERS_by_region_and_month.csv"
  ),
  row.names = FALSE
)

write.csv(
  na_results,
  file.path(
    output_csv_dir,
    "Percent_NA_by_region_and_month.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 15. COMBINE PERS, %NA AND KDE SAMPLE SIZE----
# ============================================================

analysis_monthly <- monthly_results %>%
  left_join(
    na_results,
    by = c(
      "Region",
      "Month"
    )
  ) %>%
  left_join(
    kernel_sample_sizes,
    by = c(
      "Region",
      "Month"
    )
  )


write.csv(
  analysis_monthly,
  file.path(
    output_csv_dir,
    "H2_monthly_analysis_dataset.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 16. CALCULATE REGIONAL MEDIAN PERS AND %NA----
# ============================================================

regional_summary <- analysis_monthly %>%
  group_by(
    Region
  ) %>%
  summarise(
    median_PERS = median(
      PERS,
      na.rm = TRUE
    ),
    median_percent_NA = median(
      Percent_NA,
      na.rm = TRUE
    ),
    n_months = n(),
    .groups = "drop"
  )


write.csv(
  regional_summary,
  file.path(
    output_csv_dir,
    "H2_regional_summary.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 17. ECOQO / FULMAR THRESHOLD VALUES----
# ============================================================

ecoqo_df <- data.frame(
  Region = c(
    "North Sea",
    "Faroe Plateau",
    "Iceland Shelf and Sea",
    "Canadian Eastern Arctic - West Greenland",
    "Barents Sea"
  ),
  EcoQO_FTV = c(
    51,
    40.5,
    27.6,
    14,
    22.5
  )
)


# ============================================================
# 18. FINAL HYPOTHESIS 2 DATASET----
# ============================================================

H2_analysis_df <- regional_summary %>%
  inner_join(
    ecoqo_df,
    by = "Region"
  )


write.csv(
  H2_analysis_df,
  file.path(
    output_csv_dir,
    "H2_final_analysis_dataset.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 19. KENDALL'S RANK CORRELATION----
# ============================================================

kendall_test_exact <- cor.test(
  H2_analysis_df$median_PERS,
  H2_analysis_df$EcoQO_FTV,
  method = "kendall",
  exact = TRUE
)

kendall_test_exact

# ============================================================
# 20. SAVE CORRELATION RESULTS----
# ============================================================

correlation_results <- data.frame(
  Method = "Kendall rank correlation",
  Tau = unname(
    kendall_test_exact$estimate
  ),
  P_value = kendall_test_exact$p.value,
  N_regions = nrow(
    H2_analysis_df
  )
)


write.csv(
  correlation_results,
  file.path(
    output_csv_dir,
    "H2_Kendall_correlation.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 21. FINAL DATA CHECKS----
# ============================================================

H2_analysis_df

correlation_results
