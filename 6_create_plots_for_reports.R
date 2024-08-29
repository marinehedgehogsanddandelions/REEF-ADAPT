############################################################################################
################  Create plots of variables, collate stats etc for reports  ###############
############################################################################################
# Written by G V Wood Feb 2022
# Email george.wood@flinders.edu.au for additional support

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("0_Dependencies.R")

# Load the chosen_models dataframe from the previous script
chosen_models <- read.csv("output/gdm_model_output/useful_models/chosen_models.csv", stringsAsFactors = FALSE)
#chosen_models <- chosen_models[-5,] # 5 breaks the loop, remove for now. need to sort on pval !=NA.
refid_list <- chosen_models$refid

# Assuming you have a list of directories and corresponding model names
directory_mapping <- list(
  "output/gdmmodel/connectivity_unscaled/" = "conn_unsc",
  "output/gdmmodel/connectivity_scaled/" = "conn_sc",
  "output/gdmmodel/connectivity_scaled_euc/" = "conn_sc_euc",
  "output/gdmmodel/connectivity_unscaled_euc/" = "conn_unsc_euc"
)

# Initialize a list to store the full paths of the useful models
useful_models <- list()

# Iterate over the chosen_models dataframe and create the list of useful models
for (i in 1:nrow(chosen_models)) {
  refid <- chosen_models$refid[i]
  model <- chosen_models$model[i]
  
  # Find the corresponding directory
  dir <- names(directory_mapping)[which(directory_mapping == model)]
  
  # Construct the full path to the spline file
  spline_file <- paste0(dir, "gdmspline/", "gdmspline_", refid, ".rds")
  
  # Append to the list
  useful_models <- c(useful_models, spline_file)
}

# Convert useful_models list to a vector
useful_model_splines <- unlist(useful_models)

# Print the list of useful model splines
print(useful_model_splines)

############################################################################################
################  Extract splines from GDM models on useful_models:        #################
############################################################################################

# Initialize a list to store all spline results data frames temporarily
all_spline_results <- list()

# Iterate over the useful models
for (spline_file in useful_models) {
  # Read the GDM spline file
  gdms <- readRDS(spline_file)
  
  # Extract the refid from the filename
  spline_filename <- basename(spline_file)
  refid <- gsub("gdmspline_", "", spline_filename)
  refid <- gsub(".rds", "", refid)
  
  # Check if gdms is NULL
  if (is.null(gdms)) {
    splineiresult <- data.frame("refid" = refid,
                                "x" = NA,
                                "y" = NA)
  } else {
    splineiresult <- data.frame("refid" = refid,
                                "x" = gdms$x,
                                "y" = gdms$y)
  }
  
  # Adjust column names to remove refid-specific parts
  splineiresult_colnames <- colnames(splineiresult)
  splineiresult_colnames <- sub(refid, "", splineiresult_colnames)
  splineiresult_colnames <- sub("x..", "x.", splineiresult_colnames, fixed = TRUE)
  splineiresult_colnames <- sub("y..", "y.", splineiresult_colnames, fixed = TRUE)
  colnames(splineiresult) <- splineiresult_colnames
  
  # Append the results to the list
  all_spline_results[[length(all_spline_results) + 1]] <- splineiresult
}

# Find the union of all column names
all_column_names <- unique(unlist(lapply(all_spline_results, colnames)))

# Ensure each spline result data frame has all columns
for (i in 1:length(all_spline_results)) {
  missing_columns <- setdiff(all_column_names, colnames(all_spline_results[[i]]))
  all_spline_results[[i]][missing_columns] <- NA
  all_spline_results[[i]] <- all_spline_results[[i]][, all_column_names]
}

# Combine all spline results into a single data frame
spline_results <- do.call(rbind, all_spline_results)

# Fix the refid column
spline_results$refid <- gsub(".rds", "", spline_results$refid)
spline_results$refid <- as.factor(spline_results$refid)
str(spline_results)

# Save the spline results
dir.create("output/gdm_model_output/useful_models", recursive = TRUE)
write.csv(spline_results, "output/gdm_model_output/useful_models/spline_results.csv")
saveRDS(spline_results, "output/gdm_model_output/useful_models/spline_results.rds")

#########
#### Plotting ####
#########

# Iterate over the chosen_models dataframe and create the list of useful models
useful_models <- list()

for (i in 1:nrow(chosen_models)) {
  refid <- chosen_models$refid[i]
  model <- chosen_models$model[i]
  
  # Find the corresponding directory
  dir <- names(directory_mapping)[which(directory_mapping == model)]
  
  # Construct the full path to the model file
  model_file <- paste0(dir, "gdmmodel_", refid, ".rds")
  
  # Append to the list
  useful_models <- c(useful_models, model_file)
}

# Convert useful_models list to a vector
useful_model_files <- unlist(useful_models)

# Print the list of useful model files
print(useful_model_files)

# Iterate over the useful models to create plots
for (j in 1:length(useful_model_files)) {
  skip_to_next <- FALSE
  
  tryCatch({
    refid <- refid_list[j]
    message("Processing refid: ", refid)
    
    # Derive the corresponding model file path
    model_file <- useful_model_files[j]
    
    # Read the GDM model
    modeli <- readRDS(model_file)
    
    # Extract the refid from the filename
    model_filename <- basename(model_file)
    
    # Find the corresponding directory for the model error file
    dir <- dirname(model_file)
    
    # Define the path to the model error file
    model_error_file <- paste0(dir, "/model_error/", refid, ".csv")
    
    # Extract and format the variables included in the final model
    variables_included_in_final_model <- gsub(
      pattern = paste0(refid, "\\."), 
      replacement = "", 
      x = modeli$predictors
    )
    variables_included_in_final_model <- paste(variables_included_in_final_model, collapse = ", ")
    
      modeli_error <- read.csv(model_error_file)
      colnames(modeli_error) <- gsub(paste(refid, ".", sep = ""), "", colnames(modeli_error))
      
      # Setting the y axis limits at the max amount of deviance explained for the model, so what you see is like a percentage:
      v <- c(max(modeli_error$Geographic_plusSD_Y, modeli_error$max_current_speed_min_depth_plusSD_Y, modeli_error$mean_light_at_bottom_plusSD_Y, modeli_error$BO2_tempmax_bdmin_plusSD_Y, modeli_error$BO2_temprange_bdmin_plusSD_Y, modeli_error$salinity_plusSD_Y))
      
      my_theme <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                        panel.background = element_blank(), axis.line = element_line(colour = "black"))
      
      # for geo, define the x-axis label based on the directory name
      xlab_text <- ifelse(grepl("_euc", dir), "Euclidean distance between sites", "Oceanographic connectivity")
      
      Geographic <- ggplot(modeli_error, aes(x = Geographic_fullModel_X, y = Geographic_fullModel_Y)) + 
        # Add the ribbon with a fill color
        geom_ribbon(aes(ymin = Geographic_minusSD_Y, ymax = Geographic_plusSD_Y), alpha = 0.1, fill = "green") + 
        # Add the upper and lower boundary lines with a dashed line
        geom_line(aes(y = Geographic_minusSD_Y), color = "green", linetype = "dashed") +
        geom_line(aes(y = Geographic_plusSD_Y), color = "green", linetype = "dashed") +
        # Add the main line
        geom_line(aes(y = Geographic_fullModel_Y), color = "black") + 
        my_theme + ggtitle("") + xlab(xlab_text) + ylab("Genetic turnover") + ylim(0, v)
      
      BO2_curvelmax_bdmean <- ggplot(modeli_error, aes(x = BO2_curvelmax_bdmean_fullModel_X)) + 
        geom_ribbon(aes(ymin = BO2_curvelmax_bdmean_minusSD_Y, ymax = BO2_curvelmax_bdmean_plusSD_Y), alpha = 0.1, fill = "green") + 
        geom_line(aes(y = BO2_curvelmax_bdmean_minusSD_Y), color = "green", linetype = "dashed") +
        geom_line(aes(y = BO2_curvelmax_bdmean_plusSD_Y), color = "green", linetype = "dashed") +
        geom_line(aes(y = BO2_curvelmax_bdmean_fullModel_Y), color = "black") +
        my_theme + ggtitle("") + xlab("Current speed") + ylab("Genetic turnover") + ylim(0, v)
      
      mean_light_at_bottom_full <- ggplot(modeli_error, aes(x = mean_light_at_bottom_fullModel_X)) + 
        geom_ribbon(aes(ymin = mean_light_at_bottom_minusSD_Y, ymax = mean_light_at_bottom_plusSD_Y), alpha = 0.1, fill = "green") + 
        geom_line(aes(y = mean_light_at_bottom_minusSD_Y), color = "green", linetype = "dashed") +
        geom_line(aes(y = mean_light_at_bottom_plusSD_Y), color = "green", linetype = "dashed") +
        geom_line(aes(y = mean_light_at_bottom_fullModel_Y), color = "black") +
        my_theme + ggtitle("") + xlab("Mean light") + ylab("Genetic turnover") + ylim(0, v)
      
      BO2_tempmax_bdmin <- ggplot(modeli_error, aes(x = BO2_tempmax_bdmin_fullModel_X)) + 
        geom_ribbon(aes(ymin = BO2_tempmax_bdmin_minusSD_Y, ymax = BO2_tempmax_bdmin_plusSD_Y), alpha = 0.1, fill = "green") + 
        geom_line(aes(y = BO2_tempmax_bdmin_minusSD_Y), color = "green", linetype = "dashed") +
        geom_line(aes(y = BO2_tempmax_bdmin_plusSD_Y), color = "green", linetype = "dashed") +
        geom_line(aes(y = BO2_tempmax_bdmin_fullModel_Y), color = "black") +
        my_theme + ggtitle("") + xlab("Max. temp") + ylab("Genetic turnover") + ylim(0, v)
      
      BO2_temprange_bdmin <- ggplot(modeli_error, aes(x = BO2_temprange_bdmin_fullModel_X)) + 
        geom_ribbon(aes(ymin = BO2_temprange_bdmin_minusSD_Y, ymax = BO2_temprange_bdmin_plusSD_Y), alpha = 0.1, fill = "green") + 
        geom_line(aes(y = BO2_temprange_bdmin_minusSD_Y), color = "green", linetype = "dashed") +
        geom_line(aes(y = BO2_temprange_bdmin_plusSD_Y), color = "green", linetype = "dashed") +
        geom_line(aes(y = BO2_temprange_bdmin_fullModel_Y), color = "black") +
        my_theme + ggtitle("") + xlab("Temp. range") + ylab("Genetic turnover") + ylim(0, v)
      
      BO2_salinitymax_bdmax <- ggplot(modeli_error, aes(x = BO2_salinitymax_bdmax_fullModel_X)) + 
        geom_ribbon(aes(ymin = BO2_salinitymax_bdmax_minusSD_Y, ymax = BO2_salinitymax_bdmax_plusSD_Y), alpha = 0.1, fill = "green") + 
        geom_line(aes(y = BO2_salinitymax_bdmax_minusSD_Y), color = "green", linetype = "dashed") +
        geom_line(aes(y = BO2_salinitymax_bdmax_plusSD_Y), color = "green", linetype = "dashed") +
        geom_line(aes(y = BO2_salinitymax_bdmax_fullModel_Y), color = "black") +
        my_theme + ggtitle("") + xlab("Salinity") + ylab("Genetic turnover") + ylim(0, v)
      
      Range_percentile <- ggplot(modeli_error, aes(x = Range_fullModel_X)) + 
        geom_ribbon(aes(ymin = Range_minusSD_Y, ymax = Range_plusSD_Y), alpha = 0.1, fill = "green") + 
        geom_line(aes(y = Range_minusSD_Y), color = "green", linetype = "dashed") +
        geom_line(aes(y = Range_plusSD_Y), color = "green", linetype = "dashed") +
        geom_line(aes(y = Range_fullModel_Y), color = "black") +
        my_theme + ggtitle("") + xlab("Range position") + ylab("Genetic turnover") + ylim(0, v)
      
      BO_parmean <- ggplot(modeli_error, aes(x = BO_parmean_fullModel_X)) + 
        geom_ribbon(aes(ymin = BO_parmean_minusSD_Y, ymax = BO_parmean_plusSD_Y), alpha = 0.1, fill = "green") + 
        geom_line(aes(y = BO_parmean_minusSD_Y), color = "green", linetype = "dashed") +
        geom_line(aes(y = BO_parmean_plusSD_Y), color = "green", linetype = "dashed") +
        geom_line(aes(y = BO_parmean_fullModel_Y), color = "black") +
        my_theme + ggtitle("") + xlab("PAR") + ylab("Genetic turnover") + ylim(0, v)
      
      Distance_offshore <- ggplot(modeli_error, aes(x = Distance_offshore_fullModel_X)) + 
        geom_ribbon(aes(ymin = Distance_offshore_minusSD_Y, ymax = Distance_offshore_plusSD_Y), alpha = 0.1, fill = "green") + 
        geom_line(aes(y = Distance_offshore_minusSD_Y), color = "green", linetype = "dashed") +
        geom_line(aes(y = Distance_offshore_plusSD_Y), color = "green", linetype = "dashed") +
        geom_line(aes(y = Distance_offshore_fullModel_Y), color = "black") +
        my_theme + ggtitle("") + xlab("Distance Offshore") + ylab("Genetic turnover") + ylim(0, v)
      
      Bioregion <- ggplot(modeli_error, aes(x = Bioregion_fullModel_X)) + 
        geom_ribbon(aes(ymin = Bioregion_minusSD_Y, ymax = Bioregion_plusSD_Y), alpha = 0.1, fill = "green") + 
        geom_line(aes(y = Bioregion_minusSD_Y), color = "green", linetype = "dashed") +
        geom_line(aes(y = Bioregion_plusSD_Y), color = "green", linetype = "dashed") +
        geom_line(aes(y = Bioregion_fullModel_Y), color = "black") +
        my_theme + ggtitle("") + xlab("Bioregion") + ylab("Genetic turnover") + ylim(0, v)
      
      
      # and the model fits:
      predicted <- c(modeli$predicted)
      observed <- c(modeli$observed)
      ecological <- c(modeli$ecological)
      modeli.df <- data.frame(predicted, observed, ecological)
      
      modeli_plot_pred_diss <- ggplot(modeli.df, aes(x = ecological, y = observed)) + 
        geom_point(colour = "blue") + geom_smooth(method = lm, se = FALSE, colour = "red", linewidth = 0.5) + my_theme + ggtitle("") + xlab("Predicted ecological distance") + ylab("Observed dissimilarity")
      
      modeli_plot_model_fit <- ggplot(modeli.df, aes(x = predicted, y = observed)) + 
        geom_point(colour = "blue") + geom_smooth(method = lm, se = FALSE, colour = "red", linewidth = 0.5) + my_theme + ggtitle("") + xlab("Predicted dissimilarity") + ylab("Observed dissimilarity")
      

      dir.create("output/gdm_model_output/useful_models/plots/")
      devAskNewPage(ask = FALSE)
      
      # Define the basic plots
      plots <- c("modeli_plot_pred_diss", "modeli_plot_model_fit")
      
      # Add variable-specific plots if they exist
      variable_plots <- unlist(strsplit(variables_included_in_final_model, ", "))
      variable_plots <- gsub(" ", "", variable_plots)  # Remove any extra spaces
      
      # Create the plot_objects list, checking for the existence of each plot
      plot_objects <- lapply(c(plots, variable_plots), function(x) {
        if (exists(x)) {
          return(eval(as.symbol(x)))
        } else {
          return(NULL)
        }
      })
      
      # Remove NULL entries from plot_objects
      plot_objects <- Filter(Negate(is.null), plot_objects)
      
      # Determine the number of rows required for the plot layout
      num_plots <- length(plot_objects)
      num_cols <- 2  # Fixed number of columns
      num_rows <- ceiling(num_plots / num_cols)
      
      # Construct the combined plot
      if (num_plots > 0) {
        modeli_plots <- do.call(ggarrange, c(plot_objects, list(
          labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I"), # Adjust based on the actual number of plots
          ncol = num_cols,
          nrow = num_rows
        )))
        
        # SAVE: Step 1: Call the png command to start the plot
        dir.create("output/gdm_model_output/useful_models/plots/", recursive = TRUE)
        plot_height <- num_rows * 5  # Adjust height per row (5 inches per row as an example)
        png(file = paste0("output/gdm_model_output/useful_models/plots/gdm_plots_", refid, ".png"), units = "in", width = 10, height = plot_height, res = 500) # Adjust width and height based on the number of plots
        
        plot(modeli_plots)
        
        # Step 3: Run dev.off() to create the file!
        dev.off()
      } else {
        message("No valid plots available for refid: ", refid)
      }
      
      
      ####################################
      ############# genetic structure ####
      ####################################
      
      predData1 <- readRDS(file = paste0("output/connectivity_matrix/predData/", refid, ".rds"))
      predData1 <- na.omit(predData1)
      predData <- predData1
      predData$x <- predData$NMDS1
      predData$y <- predData$NMDS2
      
      trainingData <- readRDS(file = paste0("output/connectivity_matrix/predData/trainingdata/", refid, ".rds"))
      covariates <- names(trainingData)
      covariates <- gsub(paste(refid, ".", sep = ""), "", covariates)
      
      predData <- predData %>%
        select(all_of(covariates[-1])) # do not include the site column as it does not exist in predData
      
      transRasts <- gdm.transform(model = modeli, data = predData)
      
      # Get the data from the gdm transformed rasters as a table
      rastDat <- as.data.frame(na.omit(transRasts))
      rastDat$NMDS1 <- rastDat$x
      rastDat$NMDS2 <- rastDat$y
      
      pcaSamp <- prcomp(rastDat[, 3:length(rastDat)]) # ignore first 2 columns, they have been moved to the end of the df now.
      
      rastDat$x <- predData1$x
      rastDat$y <- predData1$y
      
      # Turn the data into a raster:
      transRasts <- rasterFromXYZ(rastDat, crs = "+proj=longlat +datum=WGS84 +no_defs")
      
      # Predict the first three principle components for every cell in the rasters # note the use of the 'index' argument
      pcaRast <- predict(transRasts, pcaSamp, index = 1:3)
      
      # Scale the PCA rasters to make full use of the colour spectrum
      pcaRast[[1]] <- (pcaRast[[1]] - pcaRast[[1]]@data@min) / (pcaRast[[1]]@data@max - pcaRast[[1]]@data@min) * 255
      pcaRast[[2]] <- (pcaRast[[2]] - pcaRast[[2]]@data@min) / (pcaRast[[2]]@data@max - pcaRast[[2]]@data@min) * 255
      pcaRast[[3]] <- (pcaRast[[3]] - pcaRast[[3]]@data@min) / (pcaRast[[3]]@data@max - pcaRast[[3]]@data@min) * 255
      
      # Load the land raster
      
      land_raster <- readRDS("data/land_raster.rds")
      
      # Crop the land raster to the extent of the transformed rasters
      
      mainland_raster <- crop(land_raster, extent(transRasts))
     
      # Convert the cropped raster to polygons using the sf package
      
      mainland_poly <- as_Spatial(st_as_sf(st_as_stars(mainland_raster, as_points = FALSE)))
      # Dissolve the polygons to create a single multipart polygon
      
      mainland_poly <- st_union(st_as_sf(mainland_poly))
      sites <- read_csv("data/latlon.csv")
      my.sites <- filter(sites, REF_ID == refid)
      pts <- SpatialPointsDataFrame(coords = my.sites[, 4:3], data = my.sites,
                                    proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84 = 0,0,0"))
      # Plot the three PCA rasters simultaneously, each representing a different colour # (red, green, blue)
      
      devAskNewPage(ask = FALSE)
      png(file = paste0("output/gdm_model_output/useful_models/plots/genetic_turnover_", refid, ".png"), units = "in", width = 10, height = 7, res = 500) # The height of the plot in inches
      plotRGB(pcaRast, r = 1, g = 3, b = 2)
      plot(mainland_poly, col = "lightgrey", add = T)
      plot(pts, pch = 21, add = T, cex = 1, lwd = 1, col = "black", bg = "blue")
     # Step 3: Run dev.off() to create the file!
        
      dev.off()
    
  }, error = function(e) {
    skip_to_next <- TRUE
    message("Error processing refid: ", refid, ". Error: ", e)
  })
  if (skip_to_next) { next }
}


