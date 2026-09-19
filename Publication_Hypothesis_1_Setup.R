# ============================================================
# HYPOTHESIS 1 SETUP
# Individual 95% KDE generation for Northern fulmars
# ============================================================

# Purpose:
#   Generate individual 95% kernel density estimate (KDE)
#   rasters for each Northern fulmar tracking year during
#   the non-breeding season. These rasters are subsequently
#   used to calculate individual Plastic Exposure Risk Scores
#   (PERS).

# Analytical workflow:
#   1. Load and merge SEATRACK movement data and colony data
#   2. Filter locations to the non-breeding season (October-March)
#   3. Assign observations to tracking years
#   4. Combine Skjalfandi and Langanes as Iceland
#   5. Calculate individual KDEs using 200-km smoothing
#   6. Extract 95% KDEs
#   7. Mask land and normalise the resulting distributions
#   8. Save individual tracking-year KDE rasters

# Hypothesis 1:
#   PERS varies among Northern fulmar colonies and among
#   individuals within colonies during the non-breeding season.
# ============================================================


# ------------------------------------------------------------
# 1. Setup----
# ------------------------------------------------------------

rm(list = ls())

library(sf)
library(sp)
library(raster)
library(terra)
library(dplyr)
library(tidyverse)
library(adehabitatHR)
library(spData)
library(geosphere)


# ------------------------------------------------------------
# 2. Directories and input data----
# ------------------------------------------------------------

output_dir <- "/Users/ameydanole/Desktop/ENS_Rennes/Publication/First_hypothesis/Outputs"

input_dir <- "/Users/ameydanole/Desktop/ENS_Rennes/Rennes/argh/Amey_Danole_MS_Thesis/Ind/input_data"

locations <- readRDS(
  file.path(input_dir, "SEATRACK_FUGLA_20220307_v2.3_FA.rds")
)

summary_info <- readRDS(
  file.path(input_dir, "summaryTable.rds")
)


# ------------------------------------------------------------
# 3. Prepare movement data----
# ------------------------------------------------------------

movement_data <- merge(
  locations,
  summary_info,
  by = "ring"
)

# Retain locations from the non-breeding season (October-March)
movement_data <- movement_data %>%
  filter(
    !grepl("-04-|-05-|-06-|-07-|-08-|-09-", timestamp)
  ) %>%
  mutate(
    year = year(timestamp),
    month = month(timestamp)
  )

names(movement_data)[names(movement_data) == "ring"] <- "individ_id"


# ------------------------------------------------------------
# 4. Assign non-breeding tracking years----
# ------------------------------------------------------------

# January-June locations are assigned to the previous
# calendar year, whereas July-December locations are assigned
# to the current calendar year.

movement_data <- movement_data %>%
  mutate(
    tracking_year = ifelse(
      month < 7,
      year - 1,
      year
    )
  )


# ------------------------------------------------------------
# 5. Combine Icelandic colonies----
# ------------------------------------------------------------

movement_data <- movement_data %>%
  mutate(
    colony = ifelse(
      colony %in% c("Skjalfandi", "Langanes"),
      "Iceland",
      colony
    )
  )


# ------------------------------------------------------------
# 6. Record sample sizes----
# ------------------------------------------------------------

sample_size <- movement_data %>%
  group_by(individ_id, tracking_year) %>%
  summarise(
    Number_of_relocations = n(),
    .groups = "drop"
  )

colony_info <- movement_data %>%
  select(individ_id, colony) %>%
  distinct()

sample_size_output <- sample_size %>%
  left_join(colony_info, by = "individ_id") %>%
  select(colony, individ_id, tracking_year, Number_of_relocations)

write.csv(
  sample_size_output,
  file.path(output_dir, "Sample_size_KDE_hyp_1.csv"),
  row.names = FALSE
)


# ------------------------------------------------------------
# 7. Define land mask and WGS84 projection----
# ------------------------------------------------------------

land <- as(world, "Spatial")

proj_wgs84 <- sp::CRS(
  sp::proj4string(land)
)


# ------------------------------------------------------------
# 8. Generate individual 95% KDEs----
# ------------------------------------------------------------

for (colony_i in unique(movement_data$colony)) {
  
  colony_data <- movement_data %>%
    filter(colony == colony_i)
  
  
  for (year_i in unique(colony_data$tracking_year)) {
    
    tracks_wgs <- colony_data[
      colony_data$tracking_year == year_i,
    ]
    
    message(
      "Processing: ", colony_i, " | ", year_i,
      " | ", length(unique(tracks_wgs$individ_id)), " individuals"
    )
    
    # --------------------------------------------------------
    # 8.1 Minimum number of relocations----
    # --------------------------------------------------------
    
    if (nrow(tracks_wgs) > 4) {
      
      
      # ------------------------------------------------------
      # 8.2 Define geographic extent----
      # ------------------------------------------------------
      
      lon_min <- ifelse(
        min(tracks_wgs$lon) <= -179,
        -180,
        floor(min(tracks_wgs$lon)) - 1
      )
      
      lon_max <- ifelse(
        max(tracks_wgs$lon) >= 179,
        180,
        ceiling(max(tracks_wgs$lon)) + 1
      )
      
      lat_min <- ifelse(
        min(tracks_wgs$lat) <= -89,
        -90,
        floor(min(tracks_wgs$lat)) - 1
      )
      
      lat_max <- ifelse(
        max(tracks_wgs$lat) >= 89,
        90,
        ceiling(max(tracks_wgs$lat)) + 1
      )
      
      
      # ------------------------------------------------------
      # 8.3 Create geographic grid----
      # ------------------------------------------------------
      
      so.grid <- expand.grid(
        LON = seq(lon_min, lon_max, by = 1),
        LAT = seq(lat_min, lat_max, by = 1)
      )
      
      sp::coordinates(so.grid) <- ~LON + LAT
      sp::proj4string(so.grid) <- proj_wgs84
      
      
      # ------------------------------------------------------
      # 8.4 Define local Lambert azimuthal equal-area CRS----
      # ------------------------------------------------------
      
      mean_loc <- geosphere::geomean(
        cbind(tracks_wgs$lon, tracks_wgs$lat)
      )
      
      DgProj <- sp::CRS(
        paste0(
          "+proj=laea",
          " +lon_0=", mean_loc[1],
          " +lat_0=", mean_loc[2]
        )
      )
      
      so.grid.proj <- sp::spTransform(
        so.grid,
        CRS = DgProj
      )
      
      coords <- so.grid.proj@coords
      
      
      # ------------------------------------------------------
      # 8.5 Create 10-km KDE grid----
      # ------------------------------------------------------
      
      x_min <- min(coords[, 1]) - 1000000
      x_max <- max(coords[, 1]) + 1000000
      
      y_min <- min(coords[, 2]) - 1000000
      y_max <- max(coords[, 2]) + 1000000
      
      x_seq <- seq(x_min, x_max, by = 10000)
      y_seq <- seq(y_min, y_max, by = 10000)
      
      null.grid <- expand.grid(
        x = x_seq,
        y = y_seq
      )
      
      sp::coordinates(null.grid) <- ~x + y
      sp::gridded(null.grid) <- TRUE
      
      
      # ------------------------------------------------------
      # 8.6 Project tracking locations----
      # ------------------------------------------------------
      
      sp::coordinates(tracks_wgs) <- ~lon + lat
      sp::proj4string(tracks_wgs) <- proj_wgs84
      
      tracks <- sp::spTransform(
        tracks_wgs,
        CRS = DgProj
      )
      
      tracks$tracking_year <- factor(
        tracks@data$tracking_year
      )
      
      
      # ------------------------------------------------------
      # 8.7 Calculate individual KDEs----
      # ------------------------------------------------------
      
      # A 200-km smoothing parameter is used to account for
      # the positional uncertainty associated with GLS data.
      
      kudl <- adehabitatHR::kernelUD(
        tracks[, "individ_id"],
        grid = null.grid,
        h = 200000
      )
      
      
      # ------------------------------------------------------
      # 8.8 Convert KDEs to spatial rasters----
      # ------------------------------------------------------
      
      kern100 <- adehabitatHR::estUDm2spixdf(kudl)
      kern100_stk <- raster::stack(kern100)
      
      vud <- adehabitatHR::getvolumeUD(kudl)
      
      
      # ------------------------------------------------------
      # 8.9 Extract 95% KDEs for each individual----
      # ------------------------------------------------------
      
      individual_ids <- unique(tracks_wgs$individ_id)
      
      for (individual_j in seq_along(individual_ids)) {
        
        fud <- vud[[individual_j]]
        
        hr95 <- as.data.frame(fud)[, 1]
        hr95 <- as.numeric(hr95 <= 95)
        hr95 <- data.frame(hr95)
        
        
        # Convert 95% KDE to raster
        hr95_mod <- hr95
        
        coordinates(hr95_mod) <- coordinates(fud)
        sp::gridded(hr95_mod) <- TRUE
        
        hr95_raster <- raster(hr95_mod)
        
        
        # ----------------------------------------------------
        # 8.10 Restrict KDE to 95% home range----
        # ----------------------------------------------------
        
        rast <- kern100_stk[[individual_j]] * hr95_raster
        rast[is.na(rast)] <- 0
        
        
        # ----------------------------------------------------
        # 8.11 Normalise KDE----
        # ----------------------------------------------------
        
        rast <- rast / sum(
          raster::getValues(rast)
        )
        
        
        # ----------------------------------------------------
        # 8.12 Crop raster to non-zero extent----
        # ----------------------------------------------------
        
        rast[rast == 0] <- NA
        
        x_matrix <- is.na(as.matrix(rast))
        
        col_not_na <- which(
          colSums(x_matrix) != nrow(rast)
        )
        
        row_not_na <- which(
          rowSums(x_matrix) != ncol(rast)
        )
        
        cropped_extent <- raster::extent(
          rast,
          row_not_na[1] - 2,
          row_not_na[length(row_not_na)] + 2,
          col_not_na[1] - 2,
          col_not_na[length(col_not_na)] + 2
        )
        
        cropped <- raster::crop(
          rast,
          cropped_extent
        )
        
        cropped[is.na(rast)] <- 0
        
        
        # ----------------------------------------------------
        # 8.13 Mask land----
        # ----------------------------------------------------
        
        mask_proj <- sp::spTransform(
          land,
          DgProj
        )
        
        mask_proj_pol <- as(
          mask_proj,
          "SpatialPolygons"
        )
        
        rast_mask_na <- raster::mask(
          cropped,
          mask_proj_pol,
          inverse = TRUE
        )
        
        rast_mask <- rast_mask_na
        rast_mask[is.na(rast_mask)] <- 0
        
        
        # Renormalise after removing land
        rast_mask_sum1 <- rast_mask / sum(
          raster::getValues(rast_mask)
        )
        
        rast_mask[rast_mask == 0] <- NA
        
        rast_mask_final <- raster::mask(
          rast_mask_sum1,
          mask_proj_pol,
          inverse = TRUE
        )
        
        
        # ----------------------------------------------------
        # 8.14 Save raster in local projection----
        # ----------------------------------------------------
        
        input_raster_path <- file.path(
          output_dir,
          "input_tifs",
          paste0(
            colony_i, "_",
            individual_ids[individual_j], "_",
            year_i, ".tif"
          )
        )
        
        raster::writeRaster(
          rast_mask_final,
          filename = input_raster_path,
          format = "GTiff",
          overwrite = TRUE
        )
        
        
        # ----------------------------------------------------
        # 8.15 Save raster in WGS84----
        # ----------------------------------------------------
        
        mask_wgs84 <- raster::projectRaster(
          rast_mask_final,
          crs = proj_wgs84,
          over = FALSE
        )
        
        wgs84_raster_path <- file.path(
          output_dir,
          "unique_tifs",
          paste0(
            colony_i, "_",
            individual_ids[individual_j], "_",
            year_i, ".tif"
          )
        )
        
        raster::writeRaster(
          mask_wgs84,
          filename = wgs84_raster_path,
          format = "GTiff",
          overwrite = TRUE
        )
        
        
        # ----------------------------------------------------
        # 8.16 Save diagnostic plot----
        # ----------------------------------------------------
        
        plot_path <- file.path(
          output_dir,
          "unique_distributions",
          paste0(
            colony_i, "_",
            individual_ids[individual_j], "_",
            year_i, ".png"
          )
        )
        
        png(
          filename = plot_path,
          width = 1200,
          height = 900,
          res = 150
        )
        
        plot(mask_wgs84)
        plot(
          land,
          add = TRUE,
          col = "#66000000"
        )
        
        dev.off()
      }
    }
  }
}

# ============================================================
# 9. Calculate individual tracking-year PERS and %NA----
# ============================================================

# Individual tracking-year KDE rasters generated above are
# overlaid with the global marine floating plastic raster to
# calculate PERS and the percentage of KDE cells for which 
# the plastic distribution is undefined.

# PERS is calculated from the spatial overlap between the
# normalised KDE and the area-weighted plastic distribution.
# ============================================================


# ------------------------------------------------------------
# 9.1 Define directories and load plastic raster----
# ------------------------------------------------------------

kde_dir <- file.path(output_dir, "input_tifs")
multiplication_dir <- file.path(
  output_dir,
  "multiplication_rasters"
)

plastic_input_dir <- input_dir

plastics <- raster::raster(
  file.path(plastic_input_dir, "00_PlasticsRaster.tif")
)


# ------------------------------------------------------------
# 9.2 Normalise the plastic distribution----
# ------------------------------------------------------------

# NA values in the original plastic raster are temporarily
# treated as zero for normalisation and restored afterwards.

plastics_for_normalisation <- plastics
plastics_for_normalisation[is.na(plastics_for_normalisation)] <- 0

plastic_distribution <- plastics_for_normalisation /
  sum(raster::getValues(plastics_for_normalisation))

plastic_distribution[is.na(plastics)] <- NA


# ------------------------------------------------------------
# 9.3 Calculate latitude-dependent raster-cell areas----
# ------------------------------------------------------------

# Cell area varies with latitude. Areas are calculated in
# square metres using the WGS84 raster grid and Earth's
# mean radius.

plastic_resolution <- raster::res(plastics)
earth_radius <- 6371007.2

latitude <- raster::yFromRow(
  plastics,
  1:nrow(plastics)
)

cell_area <- (
  sin(
    pi / 180 *
      (latitude + plastic_resolution[2] / 2)
  ) -
    sin(
      pi / 180 *
        (latitude - plastic_resolution[2] / 2)
    )
) *
  (plastic_resolution[1] * pi / 180) *
  earth_radius^2

raster_cell_area <- raster::setValues(
  plastics,
  rep(cell_area, each = ncol(plastics))
)


# ------------------------------------------------------------
# 9.4 Identify individual KDE rasters----
# ------------------------------------------------------------

kde_files <- list.files(
  kde_dir,
  full.names = TRUE,
  pattern = "\\.tif$"
)


# ------------------------------------------------------------
# 9.5 Initialise output objects----
# ------------------------------------------------------------

exposure_results <- data.frame()

na_results <- data.frame()


# ------------------------------------------------------------
# 10. Calculate PERS and %NA for each tracking year----
# ------------------------------------------------------------

for (file_i in seq_along(kde_files)) {
  
  # ----------------------------------------------------------
  # 10.1 Load individual tracking-year KDE----
  # ----------------------------------------------------------
  
  kde_raster <- raster::raster(kde_files[file_i])
  
  raster_name <- kde_raster@data@names[1]
  
  
  # ----------------------------------------------------------
  # 10.2 Calculate total KDE mass before projection----
  # ----------------------------------------------------------
  
  kde_for_check <- kde_raster
  kde_for_check[is.na(kde_for_check)] <- 0
  
  kde_sum <- sum(
    raster::getValues(kde_for_check)
  )
  
  
  # ----------------------------------------------------------
  # 10.3 Project KDE onto plastic-raster grid----
  # ----------------------------------------------------------
  
  # Bilinear interpolation is used to match the KDE raster
  # to the spatial grid of the plastic distribution.
  
  kde_projected <- raster::projectRaster(
    kde_raster,
    plastics,
    method = "bilinear"
  )
  
  
  # ----------------------------------------------------------
  # 10.4 Apply latitude-dependent cell area----
  # ----------------------------------------------------------
  
  kde_area_weighted <- (
    kde_projected *
      raster_cell_area /
      100000000
  )
  
  kde_area_weighted[is.na(kde_area_weighted)] <- 0
  
  area_weighted_sum <- sum(
    raster::getValues(kde_area_weighted)
  )
  
  
  # ----------------------------------------------------------
  # 10.5 Identify KDE cells used for %NA calculation----
  # ----------------------------------------------------------
  
  kde_cells <- kde_area_weighted
  
  kde_cells[kde_cells > 0] <- 1
  kde_cells[is.na(kde_cells)] <- 0
  
  n_kde_cells <- sum(
    raster::getValues(kde_cells)
  )
  
  
  # ----------------------------------------------------------
  # 10.6 Identify undefined plastic cells----
  # ----------------------------------------------------------
  
  # Cells that were NA in the original plastic raster are
  # assigned a value of one; all defined plastic cells are
  # assigned zero.
  
  plastic_na <- plastic_distribution
  
  plastic_na[is.na(plastic_na)] <- 1
  plastic_na[plastic_na != 1] <- 0
  
  
  # ----------------------------------------------------------
  # 10.7 Calculate number of KDE cells with undefined
  # plastic distribution----
  # ----------------------------------------------------------
  
  overlapping_na_cells <- kde_cells * plastic_na
  
  n_na_cells <- sum(
    raster::getValues(overlapping_na_cells)
  )
  
  
  # ----------------------------------------------------------
  # 10.8 Calculate percentage of undefined plastic cells----
  # ----------------------------------------------------------
  
  percent_na <- (
    n_na_cells /
      n_kde_cells
  ) * 100
  
  
  # ----------------------------------------------------------
  # 10.9 Mask undefined plastic cells for PERS calculation----
  # ----------------------------------------------------------
  
  kde_area_weighted[
    is.na(plastics)
  ] <- NA
  
  
  # ----------------------------------------------------------
  # 10.10 Calculate PERS----
  # ----------------------------------------------------------
  
  exposure_overlap <- (
    kde_area_weighted *
      plastic_distribution
  )
  
  exposure_score_raster <- exposure_overlap
  exposure_score_raster[
    is.na(exposure_score_raster)
  ] <- 0
  
  exposure_score <- round(
    sum(
      raster::getValues(exposure_score_raster)
    ) * 1000000,
    4
  )
  
  
  # ----------------------------------------------------------
  # 10.11 Save exposure-overlap raster----
  # ----------------------------------------------------------
  
  output_raster <- file.path(
    multiplication_dir,
    basename(kde_files[file_i])
  )
  
  raster::writeRaster(
    kde_area_weighted,
    filename = output_raster,
    format = "GTiff",
    overwrite = TRUE
  )
  
  
  # ----------------------------------------------------------
  # 10.12 Extract colony, individual and tracking-year IDs----
  # ----------------------------------------------------------
  
  name_split <- strsplit(
    raster_name,
    "_"
  )[[1]]
  
  population <- name_split[1]
  individual <- name_split[2]
  tracking_year <- name_split[3]
  
  
  # ----------------------------------------------------------
  # 10.13 Store individual tracking-year PERS----
  # ----------------------------------------------------------
  
  exposure_results <- rbind(
    exposure_results,
    data.frame(
      population = population,
      individual = individual,
      tracking_year = tracking_year,
      kde_sum = kde_sum,
      area_weighted_sum = area_weighted_sum,
      exposure_score = exposure_score
    )
  )
  
  
  # ----------------------------------------------------------
  # 10.14 Store %NA results----
  # ----------------------------------------------------------
  
  na_results <- rbind(
    na_results,
    data.frame(
      name = raster_name,
      n_kde_cells = n_kde_cells,
      n_na_cells = n_na_cells,
      percent_na = percent_na
    )
  )
  
  message(
    "Processed ",
    file_i,
    " of ",
    length(kde_files),
    ": ",
    raster_name
  )
}


# ------------------------------------------------------------
# 11. Save individual tracking-year results----
# ------------------------------------------------------------

exposure_output_dir <- file.path(
  output_dir,
  "csv"
)

write.csv(
  exposure_results,
  file.path(
    exposure_output_dir,
    "correct_exposure_scores_by_individual_according_to_tracking_yr.csv"
  ),
  row.names = FALSE
)

write.csv(
  na_results,
  file.path(
    exposure_output_dir,
    "correct_nas_ind.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 12. Combine PERS and %NA results----
# ------------------------------------------------------------

na_results <- na_results %>%
  tidyr::separate(
    name,
    into = c("population", "individual", "tracking_year"),
    sep = "_",
    remove = FALSE
  )

merged_results <- exposure_results %>%
  left_join(
    na_results %>%
      select(
        population,
        individual,
        tracking_year,
        n_kde_cells,
        n_na_cells,
        percent_na
      ),
    by = c(
      "population",
      "individual",
      "tracking_year"
    )
  )

# ------------------------------------------------------------
# 13. Calculate within-individual variation across
#     tracking years----
# ------------------------------------------------------------

# For individuals represented in more than one tracking year,
# calculate the variance and median of tracking-year-specific
# PERS values. Median %NA is also calculated across tracking
# years.

individual_summary <- merged_results %>%
  group_by(
    population,
    individual
  ) %>%
  summarise(
    number_of_years_tracked = n(),
    intra_ind_variance_pers = var(
      exposure_score,
      na.rm = TRUE
    ),
    intra_ind_median_pers = median(
      exposure_score,
      na.rm = TRUE
    ),
    intra_ind_percent_nas = median(
      percent_na,
      na.rm = TRUE
    ),
    .groups = "drop"
  )


# ------------------------------------------------------------
# 14. Calculate colony-level H1 summary----
# ------------------------------------------------------------

# Colony-level PERS is the median of individual median PERS
# values. Between-individual variance is calculated from
# individual median PERS values, while the reported
# within-individual variance is the median of individual
# tracking-year variances.

colony_summary <- individual_summary %>%
  group_by(population) %>%
  summarise(
    pers = median(
      intra_ind_median_pers,
      na.rm = TRUE
    ),
    median_intra_ind_variance_pers = median(
      intra_ind_variance_pers,
      na.rm = TRUE
    ),
    variance_pers = var(
      intra_ind_median_pers,
      na.rm = TRUE
    ),
    percent_nas = median(
      intra_ind_percent_nas,
      na.rm = TRUE
    ),
    .groups = "drop"
  )


# ------------------------------------------------------------
# 15. Save H1 summary tables----
# ------------------------------------------------------------

write.csv(
  individual_summary,
  file.path(
    exposure_output_dir,
    "var_pers_by_ind_tracking_yr.csv"
  ),
  row.names = FALSE
)

write.csv(
  colony_summary,
  file.path(
    exposure_output_dir,
    "correct_ind_pers_by_population.csv"
  ),
  row.names = FALSE
)