
# Written by G V Wood and M. van der Mheen and M. van den Berg Feb 2022
# Email george.wood@flinders.edu.au for additional support

# set working directory and generate some objects we will need later:
wd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(wd)

libraries <- source("0_Dependencies.R")

############################################################################################
##################                   Set up Julia connectivity matrix             ##########
############################################################################################
# use script to run julia ocn (https://github.com/mheen/ocn) package in R
# Requires julia (https://julialang.org/) installed on system with JuliaConnectoR package

# path to julia OCN project
ocn_path <- paste0(wd, "/ocn/ocn/")
network_path <- paste0(wd,"/ocn/ocn/output/full_network_sparsemat.jld2")

juliaEval('using Pkg')
juliaEval(paste0('Pkg.develop(PackageSpec(path = "', ocn_path, '"))'))
juliaEval('using Ocn')
Coor <- juliaFun('Coor')
load_network <- juliaFun('load_network')
calc_fewest_edges_connectivity <- juliaFun('calc_fewest_edges_connectivity')

# load network
network <- load_network(network_path)

calculate_pairwise_distances <- function(sites) {
  n_sites <- nrow(sites)
  sources <- mapply(Coor, sites$longitude, sites$latitude)
  shortest_paths_matrix <- calc_fewest_edges_connectivity(network, sources)
  
  return(shortest_paths_matrix)
}


############################################################################################
##################                   make directories required for this script :  ##########
############################################################################################

dir.create("output/connectivity_matrix")
dir.create("output/connectivity_matrix/nMDS")
dir.create("output/connectivity_matrix/nMDS/nMDS_objects")
dir.create("output/connectivity_matrix/nMDS/nMDS_coordinates")
dir.create("output/connectivity_matrix/predData/")

############################################################################################
##################                   Generate nMDS coordinates for each species:  ##########
############################################################################################
# input raster file of range:
raster_rds_combo <- list.files(path = "output/range_percentile_rasters/transformed_percentiles_only/", pattern = ".rds", full.names = TRUE)
raster <- lapply(raster_rds_combo, readRDS)      # list of all the rasters
range_rasters <- raster::stack(raster)                # create raster stack
names(range_rasters)

latlon_spp <- read_csv("output/latlon_spp.csv")
latlon_spp_list <- split(latlon_spp, latlon_spp$refid)

metadata <- read_csv("data/metadata.csv")
relevant_studies <- subset(metadata, marker %in% c('SSR', 'SNP'))$refid # filter out studies that only have snps or msats

# Subset the list
latlon_spp_list <- latlon_spp_list[intersect(names(latlon_spp_list), relevant_studies)]
  
# Read RDS files 
provinces <- readRDS("data/MEOW/provinces.rds")
bioregions <- readRDS("data/MEOW/bioregions.rds")

# Read the shapefile using the sf package
Bioregions <- st_read("data/MEOW/meow_ecos.shp")

############################################################
#### specify to keep these objects and wipe all else after each iteration of the loop:

keep_only <- function(keep_these) {
  # Get all object names in the environment
  all_objects <- ls(envir = globalenv())
  
  # Find objects that are not in the keep_these list
  remove_these <- setdiff(all_objects, keep_these)
  
  # Remove the objects not in the keep_these list
  if (length(remove_these) > 0) {
    rm(list = remove_these, envir = globalenv())
  }
}

############################################################

# Loop through each species layer in the raster stack 
for (refid_index in 1:length(latlon_spp_list)) {
  skip_to_next <- FALSE
  tryCatch({
    # Specify objects to keep
    other_objects_needed_for_this_iteration <- c("wd", "Coor" ,"load_network", "calc_fewest_edges_connectivity", "network" , "calculate_pairwise_distances" ,"refid_index","libraries", "keep_these","keep_only","skip_to_next","raster_rds_combo", "raster", "range_rasters", "latlon_spp", "latlon_spp_list", "provinces", "bioregions", "Bioregions")
    keep_these <- c(other_objects_needed_for_this_iteration)
    # Call keep_only function
    keep_only(keep_these)
    wd
    libraries

    refid <- unique(names(latlon_spp_list)[refid_index])
    species <- unique(latlon_spp_list[[refid_index]]$spp)
    species_layer_index <- which(names(range_rasters) == species)
    Range_position <- range_rasters[[species_layer_index]]
    Range <- Range_position/Range_position # fix values to all 1 where this species is present.
    species_layer_name <- names(Range)
    
    studysites <- as.data.frame(latlon_spp_list[[refid_index]][7:8])
    pts <- SpatialPointsDataFrame(coords = studysites, data = studysites,
                                  proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
    
    rel_bioregions <- unique(raster::extract(bioregions, pts))
    
    # Subset the polygons based on the desired codes
    desired_bioregion_poly <- Bioregions[Bioregions$ECO_CODE_X %in% rel_bioregions, ]
    # Create a mask raster based on the extent of the buffered point
    Range <- crop(Range, extent(desired_bioregion_poly))
    Range <- mask(Range, desired_bioregion_poly)
    saveRDS(Range, file = paste("output/species_ranges_bioregions/", refid, ".rds", sep = ""))
    
    # Convert raster to points data frame to calculate distance
    points_df <- rasterToPoints(Range)
    points_df <- as.data.frame(points_df[, -3])
    sites <- data.frame(longitude = points_df$x, latitude = points_df$y)
    
    nmds_objects_file <- paste("output/connectivity_matrix/nMDS/nMDS_objects/", refid, ".rds", sep = "")
    nmds_coordinates_file <- paste("output/connectivity_matrix/nMDS/nMDS_coordinates/", refid, ".rds", sep = "")
    nmds_exists <- file.exists(nmds_objects_file)
    
    if (!nmds_exists) {
      
      # Call the function to calculate pairwise distances
      pairwise_distances <- calculate_pairwise_distances(sites)
      saveRDS(pairwise_distances, file = paste("output/connectivity_matrix/",refid, ".rds", sep = ""))
      
      pairwise_distances <- readRDS(paste("output/connectivity_matrix/",refid, ".rds", sep = ""))
      
      n <- nrow(pairwise_distances)
      symmetric_matrix <- matrix(NA, nrow = n, ncol = n)
      upper_tri_values <- pairwise_distances[upper.tri(pairwise_distances)]
      
      symmetric_matrix[upper.tri(symmetric_matrix)] <- upper_tri_values
      symmetric_matrix[lower.tri(symmetric_matrix)] <- t(symmetric_matrix)[lower.tri(symmetric_matrix)]
      options(scipen = 9999)
      symmetric_matrix[symmetric_matrix > 10000] <- 1000 # these either NA or are from the same cell, it is an artifact of the calculation
      reduced_symmetric_matrix <- symmetric_matrix
      reduced_symmetric_matrix[reduced_symmetric_matrix > 32] <- 33 # this is what we will call "long-distance dispersal"
      
      nMDS <- metaMDS(symmetric_matrix) # this automatically does a square root transformation
      saveRDS(nMDS, file = paste("output/connectivity_matrix/nMDS/nMDS_objects/", refid, ".rds", sep = ""))
      
      nmds_coordinates <- nMDS$points
      nmds_coordinates <- as.data.frame(nmds_coordinates)
      saveRDS(nmds_coordinates, file = paste("output/connectivity_matrix/nMDS/nMDS_coordinates/", refid, ".rds", sep = ""))
    } else {
      nmds_coordinates <- readRDS(nmds_coordinates_file)
    }
    
    ############################################################################################
    #################### Extract offshore (distance from mainland) data: #######################
    ############################################################################################
    mainland_raster <- readRDS("data/mainland_raster.rds")
    mainland_raster <- crop(mainland_raster, extent(Range))
    distance_from_mainland <- raster::distance(mainland_raster)
    extent(distance_from_mainland) <- extent(Range)
    distance_from_mainland <- mask(distance_from_mainland, Range)
    
    ### extract points for the predata:
    latlon_spp <- as.data.frame(latlon_spp)
    refidi <- refid
    latlon <- dplyr::filter(latlon_spp, refid == refidi)
    extracted_values <- data.frame(Name=latlon$Site, refid=latlon$refid, raster::extract(distance_from_mainland, pts))
    
    # Handle missing values
    na_rows <- which(is.na(extracted_values[,3]))
    for (row_index in na_rows) {
      refid <- unique(extracted_values$refid)
      site_name <- extracted_values$Name[row_index]
      site_info <- dplyr::filter(latlon, Site == site_name)
      new_cell <- move(site_info[, 3:2], distance_from_mainland)
      extracted_values$raster..extract.distance_from_mainland..pts.[row_index] <- raster::extract(distance_from_mainland, new_cell)
    }
    
    distance_offshore <- extracted_values[,-2]
    distance_offshore[2] <- distance_offshore[2]/1000
    rowname <- paste(refid, ".", "Distance_offshore", sep = "")
    colnames(distance_offshore)[2] <- c(unique(rowname))
    
    filename <- paste("output/gdmDissim/distance_offshore/", refid, ".rds", sep = "")
    saveRDS(distance_offshore, file = unique(filename))
    
    ############################################################################################
    ##################                   Make dataframe from rasters for predictions: ##########
    ############################################################################################
    environment <- load_layers(layercodes = c("BO_parmean","BO2_salinitymax_bdmax", "BO2_curvelmax_bdmean", "BO2_tempmax_bdmin", "BO2_temprange_bdmin" ), equalarea=FALSE, rasterstack=TRUE, datadir = "output/biooracle/")
    environment$Bioregion <- bioregions
    environment <- crop(environment, extent(Range))
    environment <- mask(environment, Range)
    Range_position <- crop(Range_position, extent(Range))
    environment$Range <- Range_position
    environment <- mask(environment, Range)
    environment$Distance_offshore <- distance_from_mainland/1000
    
    # extract data from the datapoints:
    points_df <- as.data.frame(rasterToPoints(environment))
    points_df$NMDS1 <- nmds_coordinates$MDS1
    points_df$NMDS2 <- nmds_coordinates$MDS2
    points_df <- points_df[, c("x", "y", "NMDS1", "NMDS2", "Range",  "Bioregion", "Distance_offshore","BO_parmean", "BO2_salinitymax_bdmax","BO2_curvelmax_bdmean", "BO2_tempmax_bdmin", "BO2_temprange_bdmin")]
    filename <- paste("output/connectivity_matrix/predData/", refid, ".rds", sep = "")
    saveRDS(points_df, file = unique(filename))
  }, error = function(e) { skip_to_next <- TRUE })
if (skip_to_next) { next }
}

