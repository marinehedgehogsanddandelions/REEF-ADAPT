# REEFADAPT APP - File preparation ####

# Written by Georgina Wood July 2023
# email: george.valentine.wood@outlook.com

###########################################
library("shiny")
library("dplyr")
library("raster")
library("gdm")
library("leaflet")
library("sf")
library("leafem")
library("mapview")
library("rgdal")


# make sure GDMs are up to date, these have just been copied and pasted from the provenance directory (gdmmodels) 
# and the entire connectivity matrix folder
# and the species specific bioregions

###########################################

# 1. Input environmental data ####
y2050_RCP_85_salinity <- raster("data/BioOracle_future_tiffs/2050_RCP85/2050AOGCM.RCP85.Benthic.Max.Depth.Salinity.Max.tif.BOv2_1.tif")
y2050_RCP_85_temp_max <- raster("data/BioOracle_future_tiffs/2050_RCP85/2050AOGCM.RCP85.Benthic.Min.Depth.Temperature.Max.tif.BOv2_1.tif")
y2050_RCP_85_temp_range <- raster("data/BioOracle_future_tiffs/2050_RCP85/2050AOGCM.RCP85.Benthic.Min.Depth.Temperature.Range.tif.BOv2_1.tif")
y2050_RCP_85_stack <- stack(y2050_RCP_85_temp_max,y2050_RCP_85_temp_range,y2050_RCP_85_salinity)
names(y2050_RCP_85_stack) <- c("BO2_tempmax_bdmin", "BO2_temprange_bdmin", "BO2_salinitymax_bdmax")

# 2. Input GDM model ####
rds_combo <- list.files(path = "output/gdmmodel/list_gdm_data/", pattern = "gdmmodel", full.names = FALSE) 
rds_combo_fullnames <- list.files( path = "output/gdmmodel/list_gdm_data/", pattern = "gdmmodel", full.names = TRUE) 
length(rds_combo)

  for(i in 1:length(rds_combo)){}
  skip_to_next <- FALSE
  tryCatch({
    modeli <- readRDS(rds_combo_fullnames[i])
    
    refnames <- rds_combo
    refidi <- gsub("gdmmodel_", "", refnames, fixed=TRUE)
    refidi <- gsub(".rds", "", refidi, fixed=TRUE)
    refidi <- refidi[i]
    
    predData1 <- readRDS(file = paste("output/connectivity_matrix/predData/", refidi, ".rds", sep = ""))
    range_cells <- predData1[,1:2]
    pts <- SpatialPointsDataFrame(coords = range_cells, data = range_cells,
                                  proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
    
    predData <- predData1
    predData$x <- predData$NMDS1
    predData$y <- predData$NMDS2
    
    trainingData <- readRDS(file = paste("output/connectivity_matrix/predData/trainingdata/", refidi, ".rds", sep = ""))
    covariates <- names(trainingData)
    covariates <- gsub(paste(refidi, ".", sep = ""), "", covariates)
    covariates
    #names(predData)
    
    predData <-  predData  %>%
      select(all_of(covariates[-1])) # do not include the site column as it does not exist in preddata
    
    
    environmenti <- predData 
    
    metadata       <- read.csv("data/metadata.csv", stringsAsFactors = TRUE) #import metadata
    
    library(dplyr)
    specific_study <- filter(metadata, refid == refidi)
    genus_name <- specific_study$Genus 
    species_name <- specific_study$Species
    species_title <- paste(genus_name,species_name,sep=" ") ## for plotting later
    full_spp <- paste(genus_name,species_name,sep="")  
    marker <- specific_study$marker
    reference <- specific_study$reference
    region <- specific_study$region
    neutral_adaptive <- specific_study$neutral.adaptive.overall.sex.linked
    number_markers <- paste("based on", specific_study$n.markers, specific_study$marker, sep = " ")

    
    # 2050 RCP 85:
    y2050_RCP_85_stacki <-  predData # this is the current data
    
    future_extracted_values <- data.frame(raster::extract(y2050_RCP_85_stack, pts)) # future data for 4 variables
    y2050_RCP_85_stacki$BO2_salinitymax_bdmax <- future_extracted_values$BO2_salinitymax_bdmax
    y2050_RCP_85_stacki$BO2_tempmax_bdmin <- future_extracted_values$BO2_tempmax_bdmin
    y2050_RCP_85_stacki$BO2_temprange_bdmin <- future_extracted_values$BO2_temprange_bdmin
    
    # make a file structure for each species:
    
    wd <- getwd() # get working directory
    mainDir <- paste(wd, "/output/transformed_gdm_rasters", sep = "")
    subDir <- full_spp
    nextdir <- paste(mainDir, "/", subDir, sep = "")
    MarkerDir <- marker
    
    ifelse(!dir.exists(file.path(mainDir, subDir)), dir.create(file.path(mainDir, subDir)), FALSE) # species directory
    ifelse(!dir.exists(file.path(nextdir, MarkerDir)), dir.create(file.path(nextdir, MarkerDir)), FALSE) # marker directory
    
    # Transform by GDM
    predData <- na.omit(predData)
    y2050_RCP_85_stacki <- na.omit(y2050_RCP_85_stacki)
    
    transRasts_current <- gdm.transform(modeli, predData)
    transRasts_2050_85 <- gdm.transform(modeli, y2050_RCP_85_stacki)

    # turn the data inrto a raster:
    test_raster <- readRDS(file = paste("output/species_ranges_bioregions/", refidi, ".rds", sep = ""))
    
    rastDat <- as.data.frame(transRasts_current)
    rastDat$NMDS1 <- rastDat$x
    rastDat$NMDS2 <- rastDat$y
    predData1 <- na.omit(predData1)
    rastDat$x <- predData1$x
    rastDat$y <- predData1$y
    transRasts_current <- rasterFromXYZ(rastDat, crs= crs(test_raster))
    
    
    rastDat <- as.data.frame(transRasts_2050_85)
    rastDat$NMDS1 <- rastDat$x
    rastDat$NMDS2 <- rastDat$y
    rastDat$x <- predData1$x
    rastDat$y <- predData1$y
    transRasts_2050_85 <- rasterFromXYZ(rastDat, crs= crs(test_raster))
    
    # save transformed rasters:
    
    dir.create("output/transformed_gdm_rasters/")
    dir.create(paste("output/transformed_gdm_rasters/", full_spp, sep = ""))
    dir.create(paste("output/transformed_gdm_rasters/", full_spp, "/", region, sep = ""))
    dir.create(paste("output/transformed_gdm_rasters/", full_spp, "/", region, "/", marker, sep = ""))
    saveRDS(transRasts_current, file = paste("output/transformed_gdm_rasters/", full_spp, "/", region,"/", marker, "/", "current","_",refidi,".rds", sep = ""))
    saveRDS(transRasts_2050_85, file = paste("output/transformed_gdm_rasters/", full_spp, "/", region,"/", marker, "/","2050_85","_",refidi, ".rds", sep = ""))
  },
  error = function(e){skip_to_next <- TRUE})
  if(skip_to_next){next}
}


#### to be updated ####
# survey gap analysis: 
## Sites with data ####
sites <- read.csv("data/R_latlon_new_studies.csv")
studyi_sites <- dplyr::filter(sites, REF_ID == refidi)
ref.coords <- unique(cbind(studyi_sites$lon, studyi_sites$lat))

ref_sites <- st_as_sf(studyi_sites, coords = c("lon", "lat"),
                      crs = crs(Range))
saveRDS (ref_sites, "output/transformed_gdm_rasters/Eckloniaradiata/SNP/Survey_sites.rds")
# extract the gdm transformed predictor values for those locations
ref.Trans.env.table <- raster::extract(transRasts,ref.coords)
ref.Trans.env.table <- ref.Trans.env.table[complete.cases(ref.Trans.env.table),]

# Calculate the similarity of each grid cell to the survey cells
# NB - This loop takes a couple of minutes
mean.similarity.region <- rep(0, length=nrow(Trans.env.table))
for(i.cell in 1:nrow(Trans.env.table))
{
  # Check if this cell has data
  if(!is.na(Trans.env.table[i.cell,ncol(Trans.env.table)]))
  {
    # loop through the reference cells, calculate similarity, add it to the tally
    for(j.cell in 1:nrow(ref.Trans.env.table))
    {
      ecol.dist.1 <- sum(abs(Trans.env.table[i.cell,c(3:ncol(Trans.env.table))] -
                               ref.Trans.env.table[j.cell,]))
      
      mean.similarity.region[i.cell] <- mean.similarity.region[i.cell] +
        (exp(-1 * (modeli$intercept + ecol.dist.1)))
    } # end for j.cells
    # Finish by dividing by the number of neighbouring cells to get the mean
    mean.similarity.region[i.cell] <- mean.similarity.region[i.cell] / nrow(ref.Trans.env.table)
  } # end if(!is.na())
} # end for i.cell
# Format the similarities into a raster an plot them
mnsim.ras <- raster(transRasts,layer=1)
mnsim.ras <- rasterize(Trans.env.table[,c(1:2)],
                       mnsim.ras,
                       field=mean.similarity.region)



# invert the similarity metric so that we are looking at a distance instead:
values( mnsim.ras) <- 1- values( mnsim.ras)
crs(mnsim.ras) = crs(Range)
extentR <- extent(Range)
mnsim.ras <- crop(mnsim.ras, extentR)

saveRDS(mnsim.ras, "output/transformed_gdm_rasters/Pocilloporadamicornis/SNP/mnsim.rds")


# plot_leaflet <- leaflet() %>%
#   addTiles() %>%
#   addCircleMarkers(data = location_plots, group = "location_plots", color = "red") %>%
#   addPopupGraphs(list(plot_list ), group = "location_plots", width = 450, height = 600)
# map style
mapstyle <- "https://api.mapbox.com/styles/v1/georgewood/cle6nzrlx000301rrklwf1173/tiles/256/{z}/{x}/{y}@2x?access_token=pk.eyJ1IjoiZ2Vvcmdld29vZCIsImEiOiJjbGNreHllNGkwZGk3M3BvNTBybG1rcmdqIn0.lhJiCvjkCp8dpXogBI5lcg"

# legend colours:
pal <- colorNumeric(c("transparent", "transparent", "transparent", "transparent",  "black"), c(0,1),
                     na.color = "transparent")

leaflet() %>%
  addTiles(urlTemplate = mapstyle, group = "Mapstyle") %>%
  addProviderTiles(providers$Esri.WorldImagery, group = "ESRI World Imagery") %>%
  addLayersControl(baseGroups = c("Mapstyle", "ESRI World Imagery")) %>%
  addRasterImage(mnsim.ras, colors = pal, opacity = 0.3) %>%
  addCircles(data = studyi_sites, group='circles',
             weight=1, radius=6000, color='red', fillColor = "red") %>%
  addMouseCoordinates() %>%
  setView(lng = 138, lat = -35, zoom = 4) %>%
  addLogo("https://static.wixstatic.com/media/998ae7_1bb005ee8e744689a521f8d787b88d78~mv2.png", position = "topright", width = 140, height = 100, offset.x = -10) %>%
  addScaleBar(position = "bottomright", options = scaleBarOptions(maxWidth = 200, metric = TRUE, imperial = FALSE))

addRasterImage(mnsim.ras, colors = pal, opacity = 0.3) %>%
  addCircles(data = studyi_sites, group='circles',
             weight=1, radius=6000, color='red', fillColor = "red")

#%>%
  addPolygons(data = rangepoly, fillColor = "white", fillOpacity = 0,
              highlightOptions = highlightOptions(color = "white", weight = 1)) 
