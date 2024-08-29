############################################################################################
############################################################################################
# Written by G V Wood and K.J. Griffin Feb 2022
# Email george.wood@flinders.edu.au for additional support

# This script uses TIF rasters to calculate range statistics for seaweeds
#  data: # Wood et al [unpub] meta-analysis metadata 
         # Species distribution predictions from Fragkopoulou, E., Serr?o, E. A., De Clerck, O., Costello, M. J., Ara?jo, M. B., Duarte, C. M., Krause-Jensen, D., & Assis, J. (2022). Global biodiversity patterns of marine forests of brown macroalgae. Global Ecology and Biogeography, 31, 636- 648. https://doi.org/10.1111/geb.13450

# set working directory:
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

source("0_Dependencies.R")
# source("1_Distributions.R")

############################################################################################
############################################################################################
# list of species names that includes only the species from the original metadata
metadata       <- read.csv("data/metadata.csv", stringsAsFactors = TRUE)
metadata$spp <- gsub(" ", "", paste(metadata$Genus, metadata$Species))
list.speciesnames <- unique(metadata$spp)

# read in species distribution rasters from Fragkopoulou, E. et al 2022.
kelps <- readRDS("output/range_percentile_rasters/rgbRaster.rds")

# read in any other taxa e.g. corals:
corals <- metadata[metadata$Class == "Cnidaria",]

# Loop through each species and read the shapefile
# Initialize a list to store the raster layers
raster_list <- list()

# Loop through each species, read the shapefile, and rasterize it
for (species in corals$spp) {
  # Construct the file path
  file_path <- paste("data/IUCN/", species, "/data_0.shp", sep = "")
  
  # Read the shapefile
  shapefile_data <- st_read(file_path)
  
  # Transform shapefile to match CRS of kelps
  shapefile_data <- st_transform(shapefile_data, CRS(projection(kelps)))
  
  # Create an empty raster with the same extent, resolution, and CRS as kelps
  empty_raster <- raster(extent(kelps), nrows=nrow(kelps), ncols=ncol(kelps), crs=crs(kelps))
  
  # Rasterize the shapefile
  raster_layer <- fasterize(shapefile_data, empty_raster)
  
  # Add the raster layer to the list
  raster_list[[species]] <- raster_layer
}

rm(shapefile_data)
rm(raster_layer)

# save coral distribution rasdters:
for (name in names(raster_list)) {
  raster <- raster_list[[name]]
  output_path <- file.path("output/range_percentile_rasters/", paste0(name, ".rds"))
  saveRDS(raster, output_path)
}

# Combine the individual raster layers into a stack (optional)
coral_raster_stack <- stack(raster_list)
rm(raster_list)

# Optionally, name the layers in the stack
names(coral_raster_stack) <- corals$spp

all_species <- stack(kelps,coral_raster_stack)
names(all_species)
relevantspecies.distributions <- subset(all_species, list.speciesnames)
rm(all_species)
rm(coral_raster_stack)
rm(kelps)

nlayers(relevantspecies.distributions) 
names.relevantspecies.distributions <- names(relevantspecies.distributions) # a list of the species in the rasterstack

# Check if any species are missing from the SDM files and add files to the rasterstack:
setdiff(list.speciesnames, names.relevantspecies.distributions) 
# "Fucusguiryi" -> plotted as fucus vesiculiosis later                                   
# "Nereialophocladia" added a shapefile based on lat/lon and depth to 15m:
Nereialophocladia <- readRDS("data/neriea_shp/Nereia.rds")
relevantspecies.distributions$Nereialophocladia <- relevantspecies.distributions$Sargassumthunbergii
values(relevantspecies.distributions$Nereialophocladia) <- NA
values(Nereialophocladia)[values(Nereialophocladia) >= -1 & values(Nereialophocladia) <= 1] <- 1
relevantspecies.distributions$Nereialophocladia <- extend(Nereialophocladia, relevantspecies.distributions$Nereialophocladia)
plot(relevantspecies.distributions$Nereialophocladia)
rm(Nereialophocladia)

calculate_percentile_range <- function(lat, mean_lat, min_lat, max_lat, is_southern) {
  if (is_southern) {
    perc_range <- ifelse(lat < mean_lat, 
                         (lat - min_lat) / (mean_lat - min_lat), 
                         (max_lat - lat) / (max_lat - mean_lat))
    lead_rear <- ifelse(lat < mean_lat, "leading", "rear")
  } else {
    perc_range <- ifelse(lat < mean_lat, 
                         (lat - min_lat) / (mean_lat - min_lat), 
                         (max_lat - lat) / (max_lat - mean_lat))
    lead_rear <- ifelse(lat < mean_lat, "rear", "leading")
  }
  return(list(perc_range = perc_range, lead_rear = lead_rear))
}


calculate_percentile_for_corals <- function(lat,median_lat, min_lat, max_lat) {
  if (lat <= median_lat) {
    perc_range <- (lat - min_lat) / (median_lat - min_lat)
  } else {
    perc_range <- (max_lat - lat) / (max_lat - median_lat)
  }
  return(perc_range)
}

# Example plotting with quadrant lines
plot_species_with_quadrants <- function(raster_data, species_name) {
  df <- as.data.frame(raster_data, xy = TRUE, na.rm = TRUE)
  colnames(df)[3] <- "layer"
  ggplot(df, aes(x = x, y = y, fill = layer)) +
    geom_raster() +
    scale_fill_viridis_c() +
    geom_vline(xintercept = -30, color = "red", linetype = "dashed") +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
    ggtitle(species_name) +
    theme_minimal()
}

# Loop for each species
for (i in names(relevantspecies.distributions)){
  sp_rasti <- relevantspecies.distributions[[i]]
  species <- names(sp_rasti)
  values(sp_rasti)[values(sp_rasti) < 1] <- NA
  cell_size <- raster::area(sp_rasti, na.rm = TRUE, weights = FALSE)
  cell_size <- cell_size[!is.na(cell_size)] # delete NAs from raster cell areas
  sum(cell_size) # this is in km2
  
  # Define range centre latitude and edges for each quadrant/hemisphere:
  sp_df <- as.data.frame(sp_rasti, xy = TRUE, na.rm = TRUE)
  
  if (species %in% corals$spp) {
    
    overall_median_lat <- median(sp_df$y, na.rm = TRUE) 
    overall_min_lat <- min(sp_df$y, na.rm = TRUE) 
    overall_max_lat <- max(sp_df$y, na.rm = TRUE) 
    sp_df$percentile_range <- mapply(calculate_percentile_for_corals, 
                                     lat = sp_df$y, 
                                     median_lat = overall_median_lat, 
                                     min_lat = overall_min_lat, 
                                     max_lat = overall_max_lat) 
    
    sp_df$lead_rear <- ifelse(sp_df$y < overall_median_lat, "rear", "lead")

  } else {
  # Assign quadrants using vectorized ifelse
  sp_df$quadrant <- ifelse(sp_df$y > 0 & sp_df$x < -30, "NHW", 
                           ifelse(sp_df$y < 0 & sp_df$x < -30, "SHW", 
                                  ifelse(sp_df$y > 0 & sp_df$x > -30, "NHE", 
                                         "SHE")))
  
  # Calculate range statistics
  quadrant_stats <- function(df, quadrant) {
    list(
      mean_lat = median(df$y[df$quadrant == quadrant], na.rm = TRUE),
      min_lat = min(df$y[df$quadrant == quadrant], na.rm = TRUE),
      max_lat = max(df$y[df$quadrant == quadrant], na.rm = TRUE)
    )
  }
  
  stats <- list(
    NHW = quadrant_stats(sp_df, "NHW"),
    NHE = quadrant_stats(sp_df, "NHE"),
    SHW = quadrant_stats(sp_df, "SHW"),
    SHE = quadrant_stats(sp_df, "SHE")
  )
  
  # Calculate percentile of each cell relative to range in each quadrant
  for (quad in unique(sp_df$quadrant)) {
    quaddf <- sp_df[sp_df$quadrant == quad, ]
    quad_stats <- stats[[quad]]
    is_southern <- quad %in% c("SHE", "SHW")
    
    percentile_info <- mapply(calculate_percentile_range, 
                              lat = quaddf$y, 
                              mean_lat = quad_stats$mean_lat, 
                              min_lat = quad_stats$min_lat, 
                              max_lat = quad_stats$max_lat, 
                              MoreArgs = list(is_southern = is_southern))
    
    sp_df$percentile_range[sp_df$quadrant == quad] <- unlist(percentile_info["perc_range", ])
    sp_df$lead_rear[sp_df$quadrant == quad] <- unlist(percentile_info["lead_rear", ])
    }
  }
  
  # Tidy up and output overall dataframe
  sp_df$species <- species
  sp_df$area <- sum(cell_size)

  # Return the range info to the raster and save
  sp_rasti[sp_rasti == 1] <- sp_df$percentile_range
  saveRDS(sp_rasti, paste("output/range_percentile_rasters/transformed_percentiles_only/", i, ".rds", sep = ""))
}

# add fucuds guiryi by duplicating fucus vesivulosis:
fucusvesiculosis <- readRDS("output/range_percentile_rasters/transformed_percentiles_only/Fucusvesiculosus.rds")
names(fucusvesiculosis) <- "Fucusguiryi"
saveRDS(fucusvesiculosis, "output/range_percentile_rasters/transformed_percentiles_only/Fucusguiryi.rds")


# Example plotting with quadrant lines
plot_species_with_quadrants <- function(raster_data, species_name) {
  df <- as.data.frame(raster_data, xy = TRUE, na.rm = TRUE)
  colnames(df)[3] <- "layer"
  ggplot(df, aes(x = x, y = y, fill = layer)) +
    geom_raster() +
    scale_fill_viridis_c() +
    geom_vline(xintercept = -30, color = "red", linetype = "dashed") +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
    ggtitle(species_name) +
    theme_minimal()
}
Macrocystispyrifera <- readRDS("output/range_percentile_rasters/transformed_percentiles_only/Macrocystispyrifera.rds")
plot_species_with_quadrants(Macrocystispyrifera, "Macrocystispyrifera")
plot_species_with_quadrants(mac, "Macrocystis pyrifera")
