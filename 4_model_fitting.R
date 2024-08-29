############################################################################################
############################################################################################

# This script runs GDM models on all studies in a loop and saves the outputs
# Written by G V Wood Feb 2022
# Email george.wood@flinders.edu.au for additional support

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("0_Dependencies.R")

############################################################################################
################  make 'sitespairtables' object needed for GDM to run correctly: ###########
############################################################################################

# THIS DATA INCLUDES THE RAW ENVIRONMENTAL VALUES EXTRACTED FROM BIO-ORACLE:
gdmDissim_list <- readRDS("output/gdmDissim/gdmDissim_list.rds") # pairwise fst data for studies 
gdmDissim_list_scaled <- readRDS("output/gdmDissim_scaled/gdmDissim_list.rds") # pairwise fst data for studies 
predData <- readRDS("output/gdmDissim/predData_uncorrelated.rds") # environmental predictor data 

metadata <- read.csv("data/metadata.csv", stringsAsFactors = TRUE)
metadata$spp <- gsub(" ", "", paste(metadata$Genus, metadata$Species))

latlon <- read.csv("data/latlon.csv")
latlon_spp <- read.csv("output/latlon_spp.csv")
latlon_spp_list <- split(latlon_spp, latlon_spp$refid)
studies <- names(latlon_spp_list) # ref ID for each study that is to be run here.

# generate site pair table objects for each study
list.gdm.data <- list()  # basic gdm table
list.gdm.data_oceandist <- list()  # gdm table with ocean distance

# generate directories needed:
dir.create("output/connectivity_matrix/predData/trainingdata")
dir.create("output/connectivity_matrix/predData/trainingdata_euclidean/")

gdmDissim_list <- gdmDissim_list[names(gdmDissim_list) %in% studies]
gdmDissim_list_scaled <- gdmDissim_list_scaled[names(gdmDissim_list_scaled) %in% studies]
predData <- predData[names(predData) %in% studies]

for (i in 1:length(studies)) {
  skip_to_next <- FALSE
  tryCatch({
    bioData1 <- as.data.frame(gdmDissim_list[i])
    bioData1_scaled <- as.data.frame(gdmDissim_list_scaled[i])
    predData1 <- as.data.frame(predData[i])
    
    names(bioData1)[1] <- "site"
    names(bioData1_scaled)[1] <- "site"
    names(predData1)[1] <- "site"
    names(predData1)[2] <- "x"
    names(predData1)[3] <- "y"
    
    # SWAP OUT THE SITE COORDS FOR THE NMDS COORDS:
    refidi <- studies[i]
    latlon <- as.data.frame(latlon_spp_list[[i]])  
    Distance_offshore <- readRDS(file = paste("output/gdmDissim/distance_offshore/", refidi, ".rds", sep = ""))
    Distance_offshore <- as.numeric(Distance_offshore[,2])
    predData1$Distance_offshore <- Distance_offshore
    
    spp <- metadata %>%
      filter(refid == refidi) %>%
      select(spp)
    
    points_df <- readRDS(file = paste("output/connectivity_matrix/predData/", refidi, ".rds", sep = ""))
    test_raster <- readRDS(file = paste("output/species_ranges_bioregions/", refidi, ".rds", sep = ""))
    quick_raster <- rasterFromXYZ(points_df, crs = crs(test_raster))
    
    library(sf)
    site_coords <- st_as_sf(predData1[,2:3], coords = c("x", "y"), crs = crs(test_raster))
    
    extracted_values <- data.frame(Name = latlon$Site, refid = latlon$refid, raster::extract(quick_raster$NMDS1, site_coords))
    
    # Handle missing values
    na_rows <- which(is.na(extracted_values[,3]))
    for (row_index in na_rows) {
      refid <- extracted_values$refid
      site_name <- extracted_values$Name[row_index]
      site_info <- filter(latlon, Site == site_name)
      new_cell <- move(site_info[, 3:2], rast(quick_raster$NMDS1))
      extracted_values$raster..extract.quick_raster.NMDS1..site_coords.[row_index] <- raster::extract(quick_raster$NMDS1, new_cell)
    }
    
    nMDS1_coords <- extracted_values[,3]
    
    extracted_values <- data.frame(Name = latlon$Site, refid = latlon$refid, raster::extract(quick_raster$NMDS2, site_coords))
    
    # Handle missing values
    na_rows <- which(is.na(extracted_values[,3]))
    for (row_index in na_rows) {
      refid <- extracted_values$refid
      site_name <- extracted_values$Name[row_index]
      site_info <- filter(latlon, Site == site_name)
      new_cell <- move(site_info[, 3:2], rast(quick_raster$NMDS2))
      extracted_values$raster..extract.quick_raster.NMDS2..site_coords.[row_index] <- raster::extract(quick_raster$NMDS2, new_cell)
    }
    
    nMDS2_coords <- extracted_values[,3]
    
    predData_euclidean_geo <- predData1
    
    predData1$x <- nMDS1_coords
    predData1$y <- nMDS2_coords
    
    saveRDS(predData1, paste("output/connectivity_matrix/predData/trainingdata/", refidi, ".rds", sep = ""))
    saveRDS(predData_euclidean_geo, paste("output/connectivity_matrix/predData/trainingdata_euclidean/", refidi, ".rds", sep = ""))
    
    # create sites pair table with the connectivity metrics as the geographic variable:
    gdmtab <- formatsitepair(bioData = bioData1, bioFormat = 3, 
                             predData = predData1, sampleSites = 1, 
                             siteColumn = "site", XColumn = "x", 
                             YColumn = "y")
    
    gdmtab <- gdmtab[!is.na(gdmtab$distance), ]  # because there were some NAs in the fst matrix which won't run in the gdm
    
    #####################################################################################
    # create sites pair table with the standard euclidean distance between cells as the geographic variable:
    gdmtab_euc <- formatsitepair(bioData = bioData1, bioFormat = 3, 
                                 predData = predData_euclidean_geo, sampleSites = 1, 
                                 siteColumn = "site", XColumn = "x", 
                                 YColumn = "y")
    
    gdmtab_euc <- gdmtab_euc[!is.na(gdmtab_euc$distance), ]  # because there were some NAs in the fst matrix which won't run in the gdm
    
    #####################################################################################
    # create sites pair table with the connectivity metrics as the geographic variable:
    gdmtab_scaled <- formatsitepair(bioData = bioData1_scaled, bioFormat = 3, 
                                    predData = predData1, sampleSites = 1, 
                                    siteColumn = "site", XColumn = "x", 
                                    YColumn = "y")
    
    gdmtab_scaled <- gdmtab_scaled[!is.na(gdmtab_scaled$distance), ]  # because there were some NAs in the fst matrix which won't run in the gdm
    
    #####################################################################################
    # create sites pair table with the standard euclidean distance between cells as the geographic variable:
    gdmtab_euc_scaled <- formatsitepair(bioData = bioData1_scaled, bioFormat = 3, 
                                        predData = predData_euclidean_geo, sampleSites = 1, 
                                        siteColumn = "site", XColumn = "x", 
                                        YColumn = "y")
    
    gdmtab_euc_scaled <- gdmtab_euc_scaled[!is.na(gdmtab_euc_scaled$distance), ]  # because there were some NAs in the fst matrix which won't run in the gdm
    
    #####################################################################################
    ################             run GDM models                          ################ 
    #####################################################################################
    
    
    # 1. Un-scaled fst with connectivity nmds values for x and y:
    tryCatch({
    rm(list = c("modeli_full_geoT", "modTest", "splinesiT"))  # remove variables that may carry over from previous model-fitting if the model fitting fails
    dir.create("output/gdmmodel/connectivity_unscaled/", recursive = TRUE, showWarnings = FALSE)
    modeli_full_geoT <- gdm(gdmtab, geo = TRUE)
    saveRDS(modeli_full_geoT, file = paste("output/gdmmodel/connectivity_unscaled/", "gdmmodel_", refidi, ".rds", sep = ""))
    
    modTest <- gdm.varImp(gdmtab, geo = TRUE, nPerm = 999, predSelect = FALSE)
    names(modTest) <- c("Model assessment", "Predictor Importance", "Predictor p-values", "Model Convergence")
    
    rownames(modTest[[2]]) <- gsub(refidi, "", rownames(modTest[[2]]))
    rownames(modTest[[3]]) <- gsub(refidi, "", rownames(modTest[[3]]))
    rownames(modTest[[4]]) <- gsub(refidi, "", rownames(modTest[[4]]))
    rownames(modTest[[2]]) <- gsub(".", "", rownames(modTest[[2]]), fixed = TRUE)
    rownames(modTest[[3]]) <- gsub(".", "", rownames(modTest[[3]]), fixed = TRUE)
    rownames(modTest[[4]]) <- gsub(".", "", rownames(modTest[[4]]), fixed = TRUE)
    
    saveRDS(modTest, file = paste("output/gdmmodel/connectivity_unscaled/", "modTest_", refidi, ".rds", sep = ""))
    
    # extract splines for custom plotting
    dir.create("output/gdmmodel/connectivity_unscaled/gdmspline/", recursive = TRUE, showWarnings = FALSE)
    splinesiT <- isplineExtract(modeli_full_geoT)
    saveRDS(splinesiT, file = paste("output/gdmmodel/connectivity_unscaled/gdmspline/", "gdmspline_", refidi, ".rds", sep = ""))
    
    
    # plot model splines with error and save csv of error values
    dir.create("output/gdmmodel/connectivity_unscaled/model_error/", recursive = TRUE, showWarnings = FALSE)
    plotUncertainty(gdmtab, sampleSites = 0.9, bsIters = 999, geo = TRUE, splines = NULL, knots = NULL, splineCol = "black", errCol = "lightgoldenrodyellow", plot.linewidth = 2.0, plot.layout = c(3, 3), parallel = FALSE, cores = 2, save = TRUE, fileName = paste("output/gdmmodel/connectivity_unscaled/model_error/", refidi, ".csv", sep = ""))
    
   
    ###################################################################################################################
    rm(list = c("modeli_full_geoT", "modTest", "splinesiT"))  # remove variables that may carry over from previous model-fitting if the model fitting fails
    }, error = function(e) { skip_to_next <- TRUE })
    
    # 2. Un-scaled fst with default euclidean distance values for geo:
    tryCatch({
    dir.create("output/gdmmodel/connectivity_unscaled_euc/", recursive = TRUE, showWarnings = FALSE)
    modeli_full_geoT <- gdm(gdmtab_euc, geo = TRUE)
    saveRDS(modeli_full_geoT, file = paste("output/gdmmodel/connectivity_unscaled_euc/", "gdmmodel_", refidi, ".rds", sep = ""))
    
    modTest <- gdm.varImp(gdmtab_euc, geo = TRUE, nPerm = 999, predSelect = FALSE)
    names(modTest) <- c("Model assessment", "Predictor Importance", "Predictor p-values", "Model Convergence")
    
    rownames(modTest[[2]]) <- gsub(refidi, "", rownames(modTest[[2]]))
    rownames(modTest[[3]]) <- gsub(refidi, "", rownames(modTest[[3]]))
    rownames(modTest[[4]]) <- gsub(refidi, "", rownames(modTest[[4]]))
    rownames(modTest[[2]]) <- gsub(".", "", rownames(modTest[[2]]), fixed = TRUE)
    rownames(modTest[[3]]) <- gsub(".", "", rownames(modTest[[3]]), fixed = TRUE)
    rownames(modTest[[4]]) <- gsub(".", "", rownames(modTest[[4]]), fixed = TRUE)
    
    saveRDS(modTest, file = paste("output/gdmmodel/connectivity_unscaled_euc/", "modTest_", refidi, ".rds", sep = ""))
    
    # extract splines for custom plotting
    dir.create("output/gdmmodel/connectivity_unscaled_euc/gdmspline/", recursive = TRUE, showWarnings = FALSE)
    splinesiT <- isplineExtract(modeli_full_geoT)
    saveRDS(splinesiT, file = paste("output/gdmmodel/connectivity_unscaled_euc/gdmspline/", "gdmspline_", refidi, ".rds", sep = ""))
    
    # plot model splines with error and save csv of error values
    dir.create("output/gdmmodel/connectivity_unscaled_euc/model_error/", recursive = TRUE, showWarnings = FALSE)
    plotUncertainty(gdmtab_euc, sampleSites = 0.9, bsIters = 999, geo = TRUE, splines = NULL, knots = NULL, splineCol = "black", errCol = "lightgoldenrodyellow", plot.linewidth = 2.0, plot.layout = c(3, 3), parallel = FALSE, cores = 2, save = TRUE, fileName = paste("output/gdmmodel/connectivity_unscaled_euc/model_error/", refidi, ".csv", sep = ""))
    
    
    ###################################################################################################################
    rm(list = c("modeli_full_geoT", "modTest", "splinesiT"))  # remove variables that may carry over from previous model-fitting if the model fitting fails
    }, error = function(e) { skip_to_next <- TRUE })
    
    # 3. Scaled fst with connectivity nmds values for x and y:
    tryCatch({
    dir.create("output/gdmmodel/connectivity_scaled/", recursive = TRUE, showWarnings = FALSE)
    modeli_full_geoT <- gdm(gdmtab_scaled, geo = TRUE)
    saveRDS(modeli_full_geoT, file = paste("output/gdmmodel/connectivity_scaled/", "gdmmodel_", refidi, ".rds", sep = ""))
    
    modTest <- gdm.varImp(gdmtab_scaled, geo = TRUE, nPerm = 999, predSelect = FALSE)
    names(modTest) <- c("Model assessment", "Predictor Importance", "Predictor p-values", "Model Convergence")
    
    rownames(modTest[[2]]) <- gsub(refidi, "", rownames(modTest[[2]]))
    rownames(modTest[[3]]) <- gsub(refidi, "", rownames(modTest[[3]]))
    rownames(modTest[[4]]) <- gsub(refidi, "", rownames(modTest[[4]]))
    rownames(modTest[[2]]) <- gsub(".", "", rownames(modTest[[2]]), fixed = TRUE)
    rownames(modTest[[3]]) <- gsub(".", "", rownames(modTest[[3]]), fixed = TRUE)
    rownames(modTest[[4]]) <- gsub(".", "", rownames(modTest[[4]]), fixed = TRUE)
    
    saveRDS(modTest, file = paste("output/gdmmodel/connectivity_scaled/", "modTest_", refidi, ".rds", sep = ""))
    
    
    # extract splines for custom plotting
    dir.create("output/gdmmodel/connectivity_scaled/gdmspline/", recursive = TRUE, showWarnings = FALSE)
    splinesiT <- isplineExtract(modeli_full_geoT)
    saveRDS(splinesiT, file = paste("output/gdmmodel/connectivity_scaled/gdmspline/", "gdmspline_", refidi, ".rds", sep = ""))
    
    
    # plot model splines with error and save csv of error values
    dir.create("output/gdmmodel/connectivity_scaled/model_error/", recursive = TRUE, showWarnings = FALSE)
    plotUncertainty(gdmtab_scaled, sampleSites = 0.9, bsIters = 999, geo = TRUE, splines = NULL, knots = NULL, splineCol = "black", errCol = "lightgoldenrodyellow", plot.linewidth = 2.0, plot.layout = c(3, 3), parallel = FALSE, cores = 2, save = TRUE, fileName = paste("output/gdmmodel/connectivity_scaled/model_error/", refidi, ".csv", sep = ""))

    
    ###################################################################################################################
    rm(list = c("modeli_full_geoT", "modTest", "splinesiT"))  # remove variables that may carry over from previous model-fitting if the model fitting fails
    }, error = function(e) { skip_to_next <- TRUE })
    
    # 4. Scaled fst with default euclidean distance values for geo:
    tryCatch({
    dir.create("output/gdmmodel/connectivity_scaled_euc/", recursive = TRUE, showWarnings = FALSE)
    modeli_full_geoT <- gdm(gdmtab_euc_scaled, geo = TRUE)
    saveRDS(modeli_full_geoT, file = paste("output/gdmmodel/connectivity_scaled_euc/", "gdmmodel_", refidi, ".rds", sep = ""))
    
    modTest <- gdm.varImp(gdmtab_euc_scaled, geo = TRUE, nPerm = 999, predSelect = FALSE)
    names(modTest) <- c("Model assessment", "Predictor Importance", "Predictor p-values", "Model Convergence")
    
    rownames(modTest[[2]]) <- gsub(refidi, "", rownames(modTest[[2]]))
    rownames(modTest[[3]]) <- gsub(refidi, "", rownames(modTest[[3]]))
    rownames(modTest[[4]]) <- gsub(refidi, "", rownames(modTest[[4]]))
    rownames(modTest[[2]]) <- gsub(".", "", rownames(modTest[[2]]), fixed = TRUE)
    rownames(modTest[[3]]) <- gsub(".", "", rownames(modTest[[3]]), fixed = TRUE)
    rownames(modTest[[4]]) <- gsub(".", "", rownames(modTest[[4]]), fixed = TRUE)
    
    saveRDS(modTest, file = paste("output/gdmmodel/connectivity_scaled_euc/", "modTest_", refidi, ".rds", sep = ""))
    
    # extract splines for custom plotting
    dir.create("output/gdmmodel/connectivity_scaled_euc/gdmspline/", recursive = TRUE, showWarnings = FALSE)
    splinesiT <- isplineExtract(modeli_full_geoT)
    saveRDS(splinesiT, file = paste("output/gdmmodel/connectivity_scaled_euc/gdmspline/", "gdmspline_", refidi, ".rds", sep = ""))
    
    # plot model splines with error and save csv of error values
    dir.create("output/gdmmodel/connectivity_scaled_euc/model_error/", recursive = TRUE, showWarnings = FALSE)
    plotUncertainty(gdmtab_euc_scaled, sampleSites = 0.9, bsIters = 999, geo = TRUE, splines = NULL, knots = NULL, splineCol = "black", errCol = "lightgoldenrodyellow", plot.linewidth = 2.0, plot.layout = c(3, 3), parallel = FALSE, cores = 2, save = TRUE, fileName = paste("output/gdmmodel/connectivity_scaled_euc/model_error/", refidi, ".csv", sep = ""))
    

    }, error = function(e) { skip_to_next <- TRUE })
  }, error = function(e) { skip_to_next <- TRUE })
  if (skip_to_next) { next }
}
