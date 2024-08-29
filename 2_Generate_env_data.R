############################################################################################
############################################################################################

# This script pulls together the environmental data and biological Fst data needed to run GDM models
# Written by G V Wood and K.J. Griffin Feb 2022
# Email george.wood@flinders.edu.au for additional support


############################################################################################
############################################################################################

# set working directory:
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("0_Dependencies.R")


############################################################################################
#######################    Load sites and environmental data:   ############################
############################################################################################

latlon <- read.csv("data/latlon.csv")             #study sites

my.sites <- data.frame(Name=c(latlon$Site), Lon=c(latlon$lon), Lat=c(latlon$lat))   # data needs to be ordered lon, lat for rasters
write.csv(my.sites, "output/my.sites.csv")

# turn your site points into a spatial points data frame:
pts <- SpatialPointsDataFrame(coords = my.sites[,c(2:3)], data = my.sites,
                              proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs"))

##### Make rasterstack of environmental covariates ####

# layers.bio2 <- list_layers( datasets="Bio-ORACLE" )       # List layers avaialble in Bio-ORACLE v2
# 
# #check coliniearity
# larger_covariate_list <- c("BO_cloudmean", "BO_calcite", "BO_damean",
#                            "BO_ph", "BO2_chlomean_bdmean", "BO_bathymean", "BO2_curvelmax_bdmean", "BO2_dissoxmean_bdmean", "BO2_ironmean_bdmean",  "BO2_phosphatemean_bdmean", "BO2_lightbotmean_bdmin", "BO2_lightbotmean_bdmax", "BO2_nitratemean_bdmean", "BO2_tempmax_bdmin", "BO2_temprange_bdmin", "BO2_salinitymax_bdmax", "BO2_tempmin_bdmin", "BO2_silicatemean_bdmean", "BO2_icecovermean_ss", "BO_parmean")
# 
# temp.present = c("BO2_curvelmax_bdmin", "BO2_lightbotmean_bdmax", "BO2_tempmax_bdmean","BO2_temprange_bdmean",
#                   "BO_salinity", "BO_parmean")
# 
# # remove redundant features:
# #calculate corrleation matrix:
# correlationMatrix <- layers_correlation(larger_covariate_list) %>% round(digits = 2)
# plot_correlation(abs(correlationMatrix))
# 
# #find attributes that are highly correlated:
# highly_correlated_pairs <- which(abs(correlationMatrix) > 0.75 & abs(correlationMatrix) < 1, arr.ind = TRUE)
# 
# 
# for (i in 1:nrow(highly_correlated_pairs)) {
#   row_index <- highly_correlated_pairs[i, 1]
#   col_index <- highly_correlated_pairs[i, 2]
# 
#   row_variable <- rownames(correlationMatrix)[row_index]
#   col_variable <- colnames(correlationMatrix)[col_index]
# 
#   cat("Variables:", row_variable, "and", col_variable, "are highly correlated.\n")
# }
# 
# variables_to_remove <- c(18,10,17) #remove phosphate, silicate and tempmin from the dataframe
# useful_variables <- larger_covariate_list[-variables_to_remove]
# 
# # Download environmental data layers
# environment <- load_layers(layercodes = useful_variables, equalarea=FALSE, rasterstack=TRUE, datadir = "output/biooracle/")
# saveRDS(environment, "output/environment.rds")

####################################################################################
#get a table describing the variables for sourcing:
# useful_variables # for some reason the subsetting using this vector won't work, so do it individually:
# df1 <- filter(layers.bio2,layers.bio2$layer_code == "BO_cloudmean")
# df2 <- filter(layers.bio2,layers.bio2$layer_code == "BO_calcite")
# df3 <- filter(layers.bio2,layers.bio2$layer_code == "BO_damean")
# df4 <- filter(layers.bio2,layers.bio2$layer_code == "BO_ph")
# df5 <- filter(layers.bio2,layers.bio2$layer_code == "BO2_chlomean_bdmean")
# df6 <- filter(layers.bio2,layers.bio2$layer_code == "BO_bathymean")
# df7 <- filter(layers.bio2,layers.bio2$layer_code == "BO2_curvelmax_bdmean")
# df8 <- filter(layers.bio2,layers.bio2$layer_code == "BO2_dissoxmean_bdmean")
# df9 <- filter(layers.bio2,layers.bio2$layer_code == "BO2_ironmean_bdmean")
# df10 <- filter(layers.bio2,layers.bio2$layer_code == "BO2_lightbotmean_bdmin")
# df11 <- filter(layers.bio2,layers.bio2$layer_code == "BO2_lightbotmean_bdmax")
# df12 <- filter(layers.bio2,layers.bio2$layer_code == "BO2_nitratemean_bdmean")
# df13 <- filter(layers.bio2,layers.bio2$layer_code == "BO2_tempmax_bdmin")
# df14 <- filter(layers.bio2,layers.bio2$layer_code == "BO2_temprange_bdmin")
# df15 <- filter(layers.bio2,layers.bio2$layer_code == "BO2_salinitymax_bdmax")
# df16 <- filter(layers.bio2,layers.bio2$layer_code == "BO2_icecovermean_ss")
# df17 <- filter(layers.bio2,layers.bio2$layer_code == "BO_parmean")
# 
# layersources <- rbind(df1, df2, df3, df4, df5, df6, df7, df8, df9, df10, df11, df12, df13, df14, df15, df16, df17)
# write.csv(layersources, file="data/environmentaldatasources.csv")

############################################################################################
#################### Extract environmental values from layers   ############################
############################################################################################

# need to ensure there are no NAs as the GDM won't run, so use points2nearestcell fn:
environment <- readRDS("output/environment.rds")

# Extract cell values at sample points
extracted_values <- data.frame(Name=my.sites$Name, raster::extract(environment, pts))
print(extracted_values)

# Find points with missing values
# Initialize an empty list to store results
na_rows_list <- list()

# Iterate over columns (excluding the first column)
for (col_name in colnames(extracted_values)[-1]) {
  na_rows <- extracted_values$Name[is.na(extracted_values[, col_name])]
  na_rows_list[[col_name]] <- na_rows
}

# Print the list of sites without data for each variable
print(na_rows_list)

# Assign nearest cell values to missing points
for (col_name in colnames(extracted_values)[-1]) {      # for each variable
  na_rows <- which(is.na(extracted_values[, col_name])) # give the row index of missing values
  
  for (row_index in na_rows) {                          # for each row index of missing values
    site_name <- extracted_values$Name[row_index]

    site_info <- my.sites[my.sites$Name == site_name, ]

    new_cells <- gecko::move(site_info[,2:3], rast(environment[[col_name]])) # move the site to the nearest location with data
    
    extracted_values[row_index, col_name] <- raster::extract(environment[[col_name]], new_cells[2,]) # specify the new value for thie missing cell

    }
}

# Print the updated dataframe
print(extracted_values)

# check no NAs:
# Find points with missing values
# Initialize an empty list to store results
na_rows_list <- list()

# Iterate over columns (excluding the first column)
for (col_name in colnames(extracted_values)[-1]) {
  na_rows <- extracted_values$Name[is.na(extracted_values[, col_name])]
  na_rows_list[[col_name]] <- na_rows
}

# Print the list of sites without data for each variable
print(na_rows_list)

saveRDS(extracted_values, "output/my.sites.environment.rds")
write.csv(extracted_values, "output/my.sites.environment.csv") 

############################################################################################
##################    Get distribution data from 001_Range_calculations.R:   ###############
############################################################################################

 #import distribution data:
raster_rds_combo <- list.files( path = "output/range_percentile_rasters/transformed_percentiles_only/", pattern = ".rds", full.names = TRUE )
raster <- lapply(raster_rds_combo, readRDS)      # list of all the rasters
percentile_stack <- stack(raster)                # create raster stack

#make df with species, refid and site gps:
metadata       <- read.csv("data/metadata.csv", stringsAsFactors = TRUE)
metadata
metadata$spp <- gsub(" ", "", paste(metadata$Genus, metadata$Species))

#merge with latlon:
 latlon <- read.csv("data/latlon.csv")
 latlon$refid <- latlon$REF_ID
 latlon_metadata <- left_join(latlon, metadata,  by = "refid")
 latlon_spp <- latlon_metadata[,2:6]
 latlon_spp$spp <- latlon_metadata$spp
 latlon_spp
 latlon_spp$Lon <- latlon_spp$lon#swap around lat and lon to lon and lat order for the extract process below
 latlon_spp$Lat <- latlon_spp$lat
 write_csv(latlon_spp, "output/latlon_spp.csv")

 # Preallocate data frame to store results
 result_data <- data.frame(Site = character(0), Value = numeric(0), stringsAsFactors = FALSE)

 # Loop through each species layer in the raster stack
 for (species_layer_index in 1:nlayers(percentile_stack)) {
   species_layer_name <- names(percentile_stack)[species_layer_index]
   
   # Find corresponding refids in metadata
  refids <- metadata %>%
     filter(spp == species_layer_name) %>%
     select(refid)
  
  refids <- paste(refids$refid, sep = "")
  
  refids <- unlist(strsplit(refids, "\\|", fixed = TRUE))
   
   # Filter my.sites to get sample sites for the current species
  sample_sites <- latlon_spp %>%
    filter(refid %in% refids)

   # Extract point data for the current species layer at the sample sites
   extracted_data <- data.frame(
     Site = sample_sites$Site,
     Value = raster::extract(percentile_stack[[species_layer_index]], sample_sites[, 7:8])
   )
   
   # Handle missing values
   na_rows <- which(is.na(extracted_data$Value))
   for (row_index in na_rows) {
     site_name <- extracted_data$Site[row_index]
     site_info <- filter(sample_sites, Site == site_name)
     new_cell <- gecko::move(site_info[, 7:8], rast(percentile_stack[[species_layer_index]]))
     extracted_data$Value[row_index] <- raster::extract(percentile_stack[[species_layer_index]], new_cell[2,])
   }
   
   # Append extracted and filled data to the result_data data frame
   result_data <- bind_rows(result_data, extracted_data)
 }
 

 # Check for duplicates in both datasets
 duplicates_result_data <- result_data[duplicated(result_data$Site) | duplicated(result_data$Site, fromLast = TRUE), ]

 saveRDS(result_data, "output/my.sites.range.rds")
 write.csv(result_data, "output/my.sites.range.csv")
 
 ############################################################################################
 #################### Extract marine bioregion data:            ############################
 ############################################################################################
 
 bioregions <- st_read(dsn = "data/MEOW/meow_ecos.shp")

 # Replace with your desired extent and resolution
 raster_template <- raster(extent(environment), resolution = res(environment))
 
 # Rasterize the shapefile using a specific attribute column, such as "bioregion_id"
 raster_bioregions <- rasterize(bioregions, raster_template, field = "ECO_CODE_X")

 saveRDS(raster_bioregions, "data/MEOW/bioregions.rds") # raster for projections
 

 raster_realms <- rasterize(bioregions, raster_template, field = "RLM_CODE")
 saveRDS(raster_realms, "data/MEOW/realms.rds") # raster for projections
 
 raster_provinces <- rasterize(bioregions, raster_template, field = "PROV_CODE")
 saveRDS(raster_provinces, "data/MEOW/provinces.rds") # raster for projections
 
 # Extract cell values at sample points
 extracted_values <- data.frame(Name=my.sites$Name, raster::extract(raster_bioregions, pts))

 # Handle missing values
 na_rows <- which(is.na(extracted_values$raster..extract.raster_bioregions..pts.)) # no NAs here.
 
 saveRDS(extracted_values, "output/my.sites.bioregions.rds")
 write_csv(extracted_values, "output/my.sites.bioregions.csv")
 
 
 ############################################################################################
 #################### merge environmental, range and bioregion data:           ##############
 ############################################################################################

 my.sites.environment  <- readRDS("output/my.sites.environment.rds")
 my.sites.range <- readRDS("output/my.sites.range.rds")
 colnames(my.sites.range)[colnames(my.sites.range) == "Site"] <- "Name" # make sure all the site names have the same name before merging data
 my.sites.bioregions <- readRDS("output/my.sites.bioregions.rds")
 
 predData <- merge(merge(my.sites.environment, my.sites.range, by = "Name"), my.sites.bioregions, by = "Name")
 colnames(predData)[colnames(predData) == "raster..extract.raster_bioregions..pts."] <- "Bioregion" # make sure all the site names have the same name before merging data
 colnames(predData)[colnames(predData) == "Value"] <- "Range" # make sure all the site names have the same name before merging data
 
 # add in lat and long columns
 predData <- merge(my.sites, predData, by = "Name")
 saveRDS(predData, "output/gdmDissim/predData_long.rds")
 
 # make it into a list
 latlon_spp
 newlatlon <- latlon_spp
 colnames(newlatlon)[colnames(newlatlon) == "Site"] <- "Name"
 predData <- merge(predData, newlatlon[, c("Name", "refid")], by.x = "Name", by.y = "Name")
 
 # put this in the correct order to help with vairation selection: - i have reduced the variables as many studies dont have many samp[le sites and models will not converge]
 # new_order <- c("Name", "Lon", "Lat", "BO2_nitratemean_bdmean", "BO2_ironmean_bdmean", "BO_ph", "BO_calcite", "BO2_chlomean_bdmean", "BO2_lightbotmean_bdmin", "BO2_lightbotmean_bdmax",  "BO_cloudmean", "BO2_icecovermean_ss",   "BO_damean" , "BO_bathymean", "BO2_dissoxmean_bdmean","BO_parmean", "Range",  "Bioregion",  "BO2_salinitymax_bdmax","BO2_tempmax_bdmin", "BO2_temprange_bdmin", "BO2_curvelmax_bdmean",  "refid")  # Replace with your variable names
 new_order <- c("Name", "Lon", "Lat", "Range",  "BO_parmean",  "BO2_salinitymax_bdmax", "Bioregion", "BO2_tempmax_bdmin", "BO2_temprange_bdmin",  "refid")  # Replace with your variable names
 
 predData  <- predData  %>%
   select(all_of(new_order))
 
 last_variable_index <- length(predData)
 predData_list =split(predData, predData$refid)
 predData_list = map(predData_list, ~ .x[,-c(last_variable_index) ]) #remove refid
 saveRDS(predData_list, "output/gdmDissim/predData.rds")
 
 # check if datapoints are correlated for each study and removes variables with no variation:

 # Initialize a list to store results
 results_list <- list()
 
 # Iterate over each element in the predData_list
 for (i in seq_along(predData_list)) {
   df <- predData_list[[i]]
   nvars <- length(df)
   
   # Extract numeric variables (excluding the first column)
   numeric_vars <- df[, 4:nvars][, sapply(df[, 4:nvars], is.numeric)]
   
   # Exclude columns with NA values
   cols_with_no_na <- colSums(is.na(numeric_vars)) == 0
   numeric_vars_no_na <- numeric_vars[, cols_with_no_na]
   
   # Exclude variables with zero standard deviation
   non_zero_sd_vars <- numeric_vars_no_na[, sapply(numeric_vars_no_na, sd) != 0]
   
   if (ncol(non_zero_sd_vars) > 1) {
     # Calculate correlations between numeric variables
     correlations <- cor(non_zero_sd_vars)
     cor_matrix_rm <- correlations  # Modify correlation matrix so only showing lower triangle
     cor_matrix_rm[upper.tri(cor_matrix_rm)] <- 0
     diag(cor_matrix_rm) <- 0
     corelationmat <- abs(cor_matrix_rm)
     
     # Remove highly correlated variables
     cols_to_keep <- !apply(corelationmat, 2, function(x) any(x > 0.8))
     reduced_variables <- non_zero_sd_vars[, cols_to_keep, drop = FALSE]
     
     new_preddata <- reduced_variables  # This has variables with no variance and high correlations removed.
     
     Name <- df[, 1]
     Lon <- df[, 2]
     Lat <- df[, 3]
     new_preddata1 <- cbind(Lat = Lat, new_preddata)
     new_preddata2 <- cbind(Lon = Lon, new_preddata1)
     new_preddata3 <- cbind(Name = Name, new_preddata2)
   } else {
     new_preddata3 <- df[, 1:3]  # If no valid variables, keep only the first three columns (Name, Lon, Lat)
   }
   
   results_list[[i]] <- new_preddata3
 }

 names(results_list) <- names(predData_list)
 
 saveRDS(results_list, "output/gdmDissim/predData_uncorrelated.rds")
 

 ############################################################################################
 ############################ Get genetic Fst data:             #############################
 ############################################################################################
 
 Fst <- read.csv("data/genetic_data.csv")
 Fst$REF_Fst <- as.numeric(Fst$REF_Fst)
 
 
 #for now on everything needs to be done in a loop for unique reference ids:
 fst_list =split(Fst, Fst$REF_ID)
 str(fst_list)
 
 # build individual Fst datasets, need to remove the initial metadata colunms and the x and y coords:
 fst_list = map(fst_list, ~ .x[,-c(1,3,4,5,6,7) ]) #remove s1 and s2 ID

 dir.create("output/fstdata")
 
 #to save these, run the code below:
 sapply(names(fst_list),
        function (x) write.table(fst_list[[x]], file=(paste("output/fstdata/", x, "_FSTdata.csv", sep=""))))
 
 
 #now reading them back in:
 temp <- list.files("output/fstdata/", pattern="*_FSTdata.csv")
 reffiles <- list.files("output/fstdata/", pattern="*_FSTdata.csv", full.names = TRUE)
 for (i in 1:length(temp)){
   assign(temp[i], read.csv(reffiles[i], sep = ""))}#now you have separate df for each study
 
############################################################################################
#########################    #make fst pairwise tables "gdmDissim":     ####################
############################################################################################

gdmDissim_list <- list()
pw_FST <- reffiles

for(i in 1:length(temp)){
  df <- assign(temp[i], read.csv(pw_FST[i], sep = ""))
  #df <- read.csv(pw_FST[i], sep = "")
  df <- df[,c(3,2,1)]
  fst_table <- df_to_pw_mat(df, from = "s1.", to = "s2.", value = "REF_Fst")
  site <- rownames(fst_table)
  gdmDissim <- cbind(site, fst_table)                                         # make Fst distance matrix "gdmDissim"
  gdmDissim_list[[length(gdmDissim_list) + 1]] <- gdmDissim                   # make a list of gdmDissim matrices

}
refids <- gsub(pattern = "_FSTdata.csv", "", temp)
names(gdmDissim_list) <- refids

saveRDS(gdmDissim_list, "output/gdmDissim/gdmDissim_list.rds")


############################################################################################
#########################    Scaling FST values to compare model fit later:     ####################
############################################################################################

# Calculate the minimum and maximum values of the column
metadata <- read.csv("data/metadata.csv")

for(i in 1:length(fst_list)){
  min_value <- min(fst_list[[i]][,1])
  max_value <- max(fst_list[[i]][,1])
  fst_list[[i]][,1] <- scale(fst_list[[i]][,1], center = min_value, scale = max_value - min_value)
  metadata[metadata$refid == fst_list[i]]$min <- min_value
  metadata[metadata$refid == fst_list[i]]$max <- max_value
}

dir.create("output/fstdata_scaled")

#to save these, run the code below:
sapply(names(fst_list),
       function (x) write.table(fst_list[[x]], file=(paste("output/fstdata_scaled/", x, "_FSTdata.csv", sep=""))))


#now reading them back in:
temp <- list.files("output/fstdata_scaled/", pattern="*_FSTdata.csv")
reffiles <- list.files("output/fstdata_scaled/", pattern="*_FSTdata.csv", full.names = TRUE)
for (i in 1:length(temp)){
  assign(temp[i], read.csv(reffiles[i], sep = ""))}#now you have separate df for each study

############################################################################################
#########################    #make fst pairwise tables "gdmDissim":     ####################
############################################################################################

gdmDissim_list <- list()
pw_FST <- reffiles

for(i in 1:length(temp)){
  df <- assign(temp[i], read.csv(pw_FST[i], sep = ""))
  #df <- read.csv(pw_FST[i], sep = "")
  df <- df[,c(3,2,1)]
  fst_table <- df_to_pw_mat(df, from = "s1.", to = "s2.", value = "REF_Fst")
  site <- rownames(fst_table)
  gdmDissim <- cbind(site, fst_table)                                         # make Fst distance matrix "gdmDissim"
  gdmDissim_list[[length(gdmDissim_list) + 1]] <- gdmDissim                   # make a list of gdmDissim matrices
  
}
refids <- gsub(pattern = "_FSTdata.csv", "", temp)
names(gdmDissim_list) <- refids
dir.create("output/gdmDissim_scaled")
saveRDS(gdmDissim_list, "output/gdmDissim_scaled/gdmDissim_list.rds")
