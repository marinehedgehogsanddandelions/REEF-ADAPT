# REEFADAPT APP ####

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
library("leaflet.extras")
library("shinyalert")
library("rgeos")
library("rmarkdown")
library("webshot")
library("shinysky")
library("htmlwidgets")
library("shinyscreenshot")

###########################################

# List variables for the sidebar:
metadata       <- read.csv("data/metadata.csv", stringsAsFactors = TRUE) #import metadata
### always leave a gap in the first row below the headings, otherwise the tool will pre-select the top refid.


genus_and_spp_list <- metadata[,8:9]
genus_and_spp_list$Region <- metadata$region
genus_and_spp_list$Marker <- metadata$marker
genus_and_spp_list$refid <- metadata$refid
unique_species_data <- genus_and_spp_list %>% distinct(Genus, Species, Region, Marker, refid)



# Define UI ####
ui <- fluidPage(
  tags$style(type = "text/css", "html, body {width:100%;height:100%}"),
  leafletOutput("ReefAdapt", height = "100vh"),
  busyIndicator(wait = 1000),
      absolutePanel(bottom = 10, left = 10, draggable = TRUE, style = "background-color: rgba(255, 255, 255, 0.7); padding: 10px; border-radius: 10px;",

      selectInput("genus_selection", "Genus:", choices = unique_species_data$Genus, selected = NULL),
      selectInput("species_selection", "Species:", choices = NULL, selected = NULL),
      selectInput("region_selection", "Region:", choices = NULL, selected = NULL),
      selectInput("marker_selection", "Marker:", choices = NULL, selected = NULL),
      selectInput("lat_selection","Latitude:",  choices = NULL, selected = NULL),
      selectInput("lon_selection","Longitude:",  choices = NULL, selected = NULL),
      selectInput("year_selection", "Year:", choices = NULL, selected = NULL),
      checkboxInput("show_advanced", "Show Advanced Options"),
      conditionalPanel(
        condition = "input.show_advanced == true",
      sliderInput("FST_threshold", "FST threshold (default is 0.05):", min = 0, max = 1, value = 0.05, step = 0.01),
      checkboxInput("show_survey_gaps", "Overlay model bounds", value = FALSE),
      checkboxInput("RCP", "RCP scenario"),
      textInput("textInputId1", "Please cite", value = "Wood, G. et al 2023. ReefAdapt: a novel tool for marine climate provenancing"),
      textInput("textInputId2", "Genetic Data Sources", value = NULL)),
      screenshotButton(label = "Take screenshot", filename = "ReefAdaptMap", server_dir = "Reports/"),
      actionButton("report", "Download Report"),
      textOutput("output_text"),)
  )
  



server <- function(input, output, session) {
  test_raster <- readRDS("output/distributions/Eckloniaradiata.rds") # so can get the crs from this as an example
  
  # Define Genus and Species outside the observeEvent blocks
  Genus <- NULL
  Species <- NULL
  Marker <- NULL
  spp <- NULL  # Define spp outside the observeEvent blocks
  clng <- NULL
  clat <- NULL
  Year <- NULL
  transRasts <- NULL
  refid <- NULL
  Region <- NULL

  
  FST_reactive <- reactiveValues(threshold = 0.05)
  
  # Update fst_threshold whenever fstThresholdInput changes
  reactive_fst_threshold <- reactive({
    FST_reactive$threshold
  })
  
  observeEvent(input$FST_threshold, {
    FST_reactive$threshold <- input$FST_threshold
    print("Updated FST threshold:")
    print(FST_reactive$threshold)
  })
  
  
  
####  
  # Define a reactiveValues to keep track of whether the click has occurred
  click_occurred <- reactiveValues(click = FALSE)
  
  # map style
  mapstyle <- "https://api.mapbox.com/styles/v1/georgewood/cle6nzrlx000301rrklwf1173/tiles/256/{z}/{x}/{y}@2x?access_token=pk.eyJ1IjoiZ2Vvcmdld29vZCIsImEiOiJjbGNreHllNGkwZGk3M3BvNTBybG1rcmdqIn0.lhJiCvjkCp8dpXogBI5lcg"
  
  # legend colours:
  pal1 <- colorNumeric(c("white", "white"), c(0,1),
                       na.color = "transparent")
  
  # Define the leaflet map only once
  leaflet_map <- leaflet() %>%
    addTiles(urlTemplate = mapstyle, group = "Mapstyle") %>%
    addProviderTiles(providers$Esri.WorldImagery, group = "ESRI World Imagery") %>%
    addLayersControl(baseGroups = c("Mapstyle", "ESRI World Imagery")) %>%
    addMouseCoordinates() %>%
    setView(lng = 138, lat = -35, zoom = 4) %>%
    addLogo("https://static.wixstatic.com/media/998ae7_1bb005ee8e744689a521f8d787b88d78~mv2.png", position = "topright", width = 140, height = 100, offset.x = -10) %>%
    addScaleBar(position = "bottomright", options = scaleBarOptions(maxWidth = 200, metric = TRUE, imperial = FALSE))
  # Render the initial leaflet map
  output$ReefAdapt <- renderLeaflet({
    leaflet_map
  })
  
  # Create reactiveValues to keep track of selected genus and species
  selected_values <- reactiveValues(genus = NULL, species = NULL)

  ## Choose species ####

  observeEvent(input$genus_selection, {
    Genus <- input$genus_selection
    spp <- paste0(Genus, Species)  # Use paste0 to concatenate without spaces
    
    # Update the species dropdown with the new choices based on the selected genus
    updateSelectInput(session, "species_selection", label = "Species:", choices = unique_species_data$Species[unique_species_data$Genus == Genus])
    
    # Check if the genus selection has changed and reset the species selection
    if (!is.null(selected_values$genus) && selected_values$genus != Genus) {
      selected_values$genus <- Genus
      selected_values$species <- NULL
      selected_values$Species <- NULL

    
      updateSelectInput(session, "species_selection", selected = NULL)
    } else {
      # If the selected genus remains the same, update the reactiveValues
      selected_values$genus <- Genus
    }
  })
  
  
  ## Choose region ####
  
  observeEvent(input$species_selection, {
    Genus <- input$genus_selection
    Species <- input$species_selection
    spp <- paste0(Genus, Species)  # Use paste0 to concatenate without spaces
    
    # Update the species dropdown with the new choices based on the selected genus
    updateSelectInput(session, "region_selection", label = "Region:", choices = unique_species_data$Region[unique_species_data$Species == Species & unique_species_data$Genus == Genus])
    
    # Check if the genus selection has changed and reset the species selection
    if (!is.null(selected_values$species) && selected_values$species != Species) {
      selected_values$species <- Species
      selected_values$species <- NULL
      selected_values$Species <- NULL
      updateSelectInput(session, "region_selection", selected = NULL)
    } else {
      # If the selected genus remains the same, update the reactiveValues
      selected_values$species <- Species
    }
  })
  ## Choose marker ####
  observeEvent(input$region_selection, {
    
    # Reset all other inputs when species selection changes
    updateSelectInput(session, "marker_selection", label = "Marker:", choices = NULL, selected = NULL)
    updateSelectInput(session, "lat_selection", label = "Latitude:", choices = NULL, selected = NULL)
    updateSelectInput(session, "lon_selection", label = "Longitude:", choices = NULL, selected = NULL)
    updateSelectInput(session, "year_selection", label = "Year:", choices = NULL, selected = NULL)

    Region <- input$region_selection
    spp <- paste0(Genus, Species)  # Use paste0 to concatenate without spaces

    # Update the marker dropdown with the new choices based on the selected species
    updateSelectInput(session, "marker_selection", label = "Marker:", choices = unique_species_data$Marker[unique_species_data$Region == Region], selected = NULL)
  })
  
  
  # 1.  Input distribtuion ####
  observeEvent(c(input$marker_selection), {
    req(input$species_selection, input$genus_selection, input$region_selection, input$marker_selection)  # Wait for all selections to be made
    Genus <- input$genus_selection
    Species <- input$species_selection
    Marker <- input$marker_selection
    Region <- input$region_selection
    spp <- paste0(Genus, Species)  # Use paste0 to concatenate without spaces
    
    study_path <- paste("output","/", "transformed_gdm_rasters", "/",spp, "/", Region, "/", Marker, "/", sep = ""  )
    print(study_path)
    rds_combo <- list.files(path = study_path, pattern = "current_", full.names = FALSE) 
    refid <- gsub("current_", "", rds_combo)
    refid <- gsub(".rds", "", refid)
    
    
    # 1.  Map details ####
    # species range 
    raster_file <- paste("output/species_ranges_bioregions/", refid, ".rds", sep = "") # Modify the filename format as per your needs
    print(raster_file)  # Print raster_file to check the constructed filename
    Range <- readRDS(raster_file)
    Range <- Range/Range
    rangepoly <- rasterToPolygons(Range)
    rangepoly <- gUnaryUnion(rangepoly)

    ### add to map ####
    
    leafletProxy('ReefAdapt') %>%
      clearImages() %>%
       clearPopups() 

    leafletProxy('ReefAdapt') %>%
       addPolygons(data = rangepoly, fillColor = "white", fillOpacity = 0.3,
                  highlightOptions = highlightOptions(color = "white", weight = 1)) 


    # #### survey gaps ####
    observeEvent(input$show_survey_gaps, {
      if (input$show_survey_gaps) {
    Marker <- input$marker_selection
    Genus <- input$genus_selection
    Species <- input$species_selection
    spp <- paste0(Genus, Species)  # Use paste0 to concatenate without spaces
  
    filepath1 <- paste("output/transformed_gdm_rasters/", spp, "/", Region, "/", Marker, sep = "")
    SS1 <- list.files(path = (filepath1), pattern = "Survey_sites")
    SS2 <- paste(filepath1, "/", SS1, sep = "")
    ref_sites <- readRDS(SS2)

    gap1 <- list.files(path = (filepath1), pattern = "mnsim")
    gap2 <- paste(filepath1, "/", gap1, sep = "")
    gaps <- readRDS(gap2)
    extentR <- extent(Range)
    gaps <- crop(gaps, extentR) # must be the same extent or it wont show up

    pal <- colorNumeric(c("white","black"), values(gaps),
                        na.color = "transparent")
    leafletProxy('ReefAdapt') %>%
      clearShapes() %>%
      addRasterImage(gaps, colors = pal, opacity = 0.3) %>%
      addCircles(data = ref_sites, group='circles',
                 weight=1, radius=6000, color='red', fillColor = "red")%>%
      addPolygons(data = rangepoly, fillColor = "white", fillOpacity = 0,
                  highlightOptions = highlightOptions(color = "white", weight = 1)) 

      }
   })

    
   # 2. Observe mouse clicks ####
    observeEvent(input$ReefAdapt_click, {

      click <- input$ReefAdapt_click
      clat <- click$lat
      clng <- click$lng
      

      leafletProxy('ReefAdapt') %>%
        clearGroup('circles') %>%
        clearImages() %>%
        addCircles(lng=clng, lat=clat, group='circles',
                   weight=1, radius=6000, color='black', fillColor='green',
                   fillOpacity=0.2, opacity=1) 
      
      # Set click_occurred to TRUE when the click event occurs
      click_occurred$click <- TRUE
      
      # After the click event, update the choices in the "year" dropdown
      updateSelectInput(session, "year_selection", label = "Select a year:",
                        choices = c("", "2020", "2050"))
      updateSelectInput(session, "lat_selection", label = "Latitude",
                        choices = clat)
      updateSelectInput(session, "lon_selection", label = "Longitude",
                        choices = clng)
      
      # Print the coordinates to the console
      cat("Clicked coordinates: ", clng, ", ", clat, "\n")
      
      # Retrieve the clicked coordinates from the reactive object
      
      rest.site <- data.frame(lon = clng, lat = clat)
      rest.site <- SpatialPointsDataFrame(coords = rest.site, proj4string = crs(Range), data = rest.site)
      
      print(rest.site)
      
    })
    
    
    # Observe the year dropdown
    observe({
      # If the click has not occurred yet, do not update the choices
      if (!click_occurred$click) {
        updateSelectInput(session, "year_selection", label = "Year:",
                          choices = NULL)
      }
    })
    
    # 4. Similarity analysis ####
    ## Choose year ####
    observeEvent(input$year_selection, {
      year <- input$year_selection
      click <- input$ReefAdapt_click
      Marker <- input$marker_selection
      Genus <- input$genus_selection
      Species <- input$species_selection
      spp <- paste0(Genus, Species)  # Use paste0 to concatenate without spaces
      #rest.site <- data.frame(lon = as.numeric(146.56649), lat = as.numeric(-39.06199))
      rest.site <- data.frame(lon = as.numeric(input$lon_selection), lat = as.numeric(input$lat_selection))
      str(rest.site)
      
      if (!anyNA(rest.site)){
        rest.site1 <- SpatialPointsDataFrame(coords = rest.site, proj4string = crs(test_raster), data = rest.site)
        
        
        # Observe the year dropdown
        observe({
          # If the click has not occurred yet, do not update the choices
          if (year == "2020") {
            filepath1 <- paste("output/transformed_gdm_rasters/", spp, "/", Region, "/", Marker, sep = "")
            transRasts1 <- list.files(path = (filepath1), pattern = "current")
            filepath2 <- paste(filepath1, "/", transRasts1, sep = "")
            transRasts <- readRDS(filepath2)

            
            # Extract the transformed environmental values for these focal locations
            # check both in same projection
            rest.site1 <- spTransform(rest.site1, crs(transRasts))
            focal.trans <- raster::extract(transRasts, rest.site1)
            
            # check this is running
            print(focal.trans)

            
            source <- metadata %>%
              filter(refid == refid)
            source <- paste(source$reference)
            updateTextInput(session,"textInputId2", value = source)
            
            modelipath <- paste("output/gdmmodel/list_gdm_data/", "gdmmodel_",refid, ".rds", sep = "")
            modeli <- readRDS(modelipath)
            
            # put the values from the transformed layers in a table for easy analysis
            Trans.env.table <- as.matrix(transRasts)
            print(Trans.env.table)
            col.longs<-xFromCol(transRasts)
            row.lats<-yFromRow(transRasts)
            Cell_Long<-rep(col.longs, times=nrow(transRasts))
            Cell_Lat<-rep(row.lats, each=ncol(transRasts), times=1)
            Trans.env.table<-cbind(Cell_Long, Cell_Lat, Trans.env.table)
            Trans.env.table <- Trans.env.table[complete.cases(Trans.env.table),]
            
            # now calculate the similarity of all other grid cells to each of these focal locations
            similarity.focal.pt.1 <- rep(0, length=nrow(Trans.env.table))

            for(i.cell in 1:nrow(Trans.env.table))
            {
              ecol.dist.1 <- sum(abs(Trans.env.table[i.cell,c(3:ncol(Trans.env.table))] - focal.trans[1,]))
              similarity.focal.pt.1[i.cell] <- exp(-1 * (modeli$intercept + ecol.dist.1))
        
                } # end for i.cell
            
            # Format the similarities into a raster an plot them
            focal.pt.1.ras <- raster(transRasts,layer=1)
            focal.pt.1.ras <- rasterize(Trans.env.table[,c(1:2)], focal.pt.1.ras,
                                        field=similarity.focal.pt.1)
            focal.pt.1.ras_FST_predicted <- focal.pt.1.ras
            
            
          #  if (refid != "ref_225") { # dont minus it from 1 as have already unscaled the data, for this refid only
            values(focal.pt.1.ras_FST_predicted) <- 1-values(focal.pt.1.ras)
           # }
            
            
            
            #### un-scaled scaled data (at this moment this is just for a. kenti so directly coding it below. if adding more, will need a spreadsheet with refid, max and min values.)
            if (refid == "ref_225") {
              max_value <- 0.042433
              min_value <- 0.008971
              values(focal.pt.1.ras_FST_predicted) <-  (values(focal.pt.1.ras_FST_predicted) * (max_value - min_value)) + min_value
            }
            
            
            if (refid == "ref_226") {
              max_value <- 0.17369176
              min_value <- 0.03351753
              values(focal.pt.1.ras_FST_predicted) <-  (values(focal.pt.1.ras_FST_predicted) * (max_value - min_value)) + min_value
            }
            
            if (refid == "ref_176") {
              max_value <- 0.148
              min_value <- 0.005
              values(focal.pt.1.ras_FST_predicted) <-  (values(focal.pt.1.ras_FST_predicted) * (max_value - min_value)) + min_value
            }
            
            ########
            
            
            map_values1 <- focal.pt.1.ras_FST_predicted
            
            r1NaM <- is.na(as.matrix(map_values1))
            colNotNA <- which(colSums(r1NaM) != nrow(map_values1))
            rowNotNA <- which(rowSums(r1NaM) != ncol(map_values1))
            r3Extent <- extent(map_values1, rowNotNA[1], rowNotNA[length(rowNotNA)],
                               colNotNA[1], colNotNA[length(colNotNA)])
            map_values <- crop(map_values1, r3Extent)
            print(!is.na(map_values))
            

        
            # Create a logical condition for subsetting
            upper_limit <- as.numeric(FST_reactive$threshold)
           
            # round numbers to the closest round value:
            map_values <- round(map_values, digits = 2)
            
            values(map_values)[values(map_values) <= upper_limit] <- upper_limit
            values(map_values)[values(map_values) > upper_limit] <- NA
           
            #in case there were no cells within that threshold, transform the target cell so at least it shows up as reasonable.
            # note that this is not appropriate for future cases and will only be used for current scenarios
            cell_indices <- cellFromXY(map_values, rest.site1)
            map_values[cell_indices] <- upper_limit # sets the target square to upper limit so it still shows up. could make this colour coded so you know it is actually above the threshold.
          
            
            map_values2 <- is.na(as.matrix(map_values))
            colNotNA <- which(colSums(map_values2)  != nrow(map_values))
            rowNotNA <- which(rowSums(map_values2) != ncol(map_values))
            map_values3Extent <- extent(map_values, rowNotNA[1], rowNotNA[length(rowNotNA)],
                               colNotNA[1], colNotNA[length(colNotNA)])
            map_values3 <- crop(map_values, map_values3Extent)

            cropped_poly <- rasterToPolygons(map_values3)
            cropped_poly <- gUnaryUnion(cropped_poly)
            # Simplify the merged polygon if needed
            simplified_polygon <- gSimplify(cropped_poly, tol = 0.001)

            
            leafletProxy('ReefAdapt') %>%
              clearImages() %>%
              addPolygons(data =  cropped_poly, fillColor = "limegreen", fillOpacity = 0.3,
                          highlightOptions = highlightOptions(color = "limegreen", weight = 1)
              )

            
            # 
            # breaks <- c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)
            # colours <- c("green", "transparent", "transparent", "transparent", "transparent", "transparent", "transparent", "transparent", "transparent", "transparent")
            # pal <- colorBin(colours, bins = breaks, na.color = "transparent")
            # 
            # leafletProxy('ReefAdapt') %>%
            #   clearImages()
            # 
            # leafletProxy('ReefAdapt') %>%
            #   addRasterImage(map_values, colors = pal, opacity = 0.4)
            # print(summary(map_values))
          }
          
          observe({
            # If the click has not occurred yet, do not update the choices
            if (year == "2050") {
              
              filepath1 <- paste("output/transformed_gdm_rasters/", spp, "/", Region, "/", Marker, sep = "") # you need this as you are comparing the site in the future belw to normal conditions today
              transRasts1 <- list.files(path = (filepath1), pattern = "current")
              filepath2 <- paste(filepath1, "/", transRasts1, sep = "")
              transRasts <- readRDS(filepath2)
              
              
              filepath1 <- paste("output/transformed_gdm_rasters/", spp, "/", Region, "/", Marker, sep = "")
              transRasts1 <- list.files(path = (filepath1), pattern = "2050_85_")
              filepath2 <- paste(filepath1, "/", transRasts1, sep = "")
              transRasts2050 <- readRDS(filepath2)
              
              rest.site1_2050 <- spTransform(rest.site1, crs(transRasts2050))
              focal.trans_2050 <- raster::extract(transRasts2050, rest.site1_2050)

              varname <- names(transRasts[[3]])
              source <- metadata %>%
                filter(refid == refid)
              source <- paste(source$reference)
              updateTextInput(session,"textInputId2", value = source)
              
              
              modelipath <- paste("output/gdmmodel/list_gdm_data/", "gdmmodel_",refid, ".rds", sep = "")
              modeli <- readRDS(modelipath)
              
              # put the values from the transformed layers in a table for easy analysis
              Trans.env.table <- as.matrix(transRasts)
              col.longs<-xFromCol(transRasts)
              row.lats<-yFromRow(transRasts)
              Cell_Long<-rep(col.longs, times=nrow(transRasts))
              Cell_Lat<-rep(row.lats, each=ncol(transRasts), times=1)
              Trans.env.table<-cbind(Cell_Long, Cell_Lat, Trans.env.table)
              Trans.env.table <- Trans.env.table[complete.cases(Trans.env.table),]
              
              # now calculate the similarity of all other grid cells to each of these focal locations
              similarity.focal.pt.1 <- rep(0, length=nrow(Trans.env.table))
 
              for(i.cell in 1:nrow(Trans.env.table))
              {
                ecol.dist.1 <- sum(abs(Trans.env.table[i.cell,c(3:ncol(Trans.env.table))] - focal.trans_2050[1,]))
                similarity.focal.pt.1[i.cell] <- exp(-1 * (modeli$intercept + ecol.dist.1))
              
                
              } # end for i.cell
              
              # Format the similarities into a raster an plot them
              focal.pt.1.ras_2050 <- raster(transRasts,layer=1)
              focal.pt.1.ras_2050 <- rasterize(Trans.env.table[,c(1:2)], focal.pt.1.ras_2050,
                                          field=similarity.focal.pt.1)
              focal.pt.1.ras_FST_predicted_2050 <- focal.pt.1.ras_2050
              
              #  if (refid != "ref_225") {
              values(focal.pt.1.ras_FST_predicted_2050) <- 1-values(focal.pt.1.ras_2050)
            #  }
              #### un-scaled scaled data (at this moment this is just for a. kenti so directly coding it below. if adding more, will need a spreadsheet with refid, max and min values.)
              if (refid == "ref_225") {
                max_value <- 0.042433
                min_value <- 0.008971
                values(focal.pt.1.ras_FST_predicted_2050) <-  (values(focal.pt.1.ras_FST_predicted_2050) * (max_value - min_value)) + min_value
              }
              
              map_values1_2050 <- focal.pt.1.ras_FST_predicted_2050
              r1NaM <- is.na(as.matrix(map_values1_2050))
              colNotNA <- which(colSums(r1NaM) != nrow(map_values1_2050))
              rowNotNA <- which(rowSums(r1NaM) != ncol(map_values1_2050))
              r3Extent <- extent(map_values1_2050, rowNotNA[1], rowNotNA[length(rowNotNA)],
                                 colNotNA[1], colNotNA[length(colNotNA)])
              map_values <- crop(map_values1_2050, r3Extent)
            
              # Create a logical condition for subsetting
              upper_limit <- as.numeric(FST_reactive$threshold)
              
              # round numbers to the closest round value:
              map_values1_2050 <- round(map_values1_2050, digits = 2)
              
              values(map_values1_2050)[values(map_values1_2050) <= upper_limit] <- upper_limit
              values(map_values1_2050)[values(map_values1_2050) > upper_limit] <- NA
              
              values_string <- values(map_values1_2050)
              
              if (all(is.na(values_string))) {
                shinyalert::shinyalert(
                  title = "Alert!",
                  text = "We cannot find populations adapted to these conditions at this time.

                  Please note that ReefAdapt is currently in beta testing mode",
                  type = "info"
                )
                
              } else {
              
              
              map_values2 <- is.na(as.matrix(map_values1_2050))
              colNotNA <- which(colSums(map_values2)  != nrow(map_values))
              rowNotNA <- which(rowSums(map_values2) != ncol(map_values))
              map_values3Extent <- extent(map_values, rowNotNA[1], rowNotNA[length(rowNotNA)],
                                          colNotNA[1], colNotNA[length(colNotNA)])
              map_values3 <- crop(map_values1_2050, map_values3Extent)

                cropped_poly <- rasterToPolygons(map_values3)
                cropped_poly <- gUnaryUnion(cropped_poly)
                # Simplify the merged polygon if needed
                simplified_polygon <- gSimplify(cropped_poly, tol = 0.001)

                leafletProxy('ReefAdapt') %>%
                  clearImages()

                leafletProxy('ReefAdapt') %>%
                  addPolygons(data =  cropped_poly, fillColor = "red", fillOpacity = 0.3,
                              highlightOptions = highlightOptions(color = "red", weight = 1)
                  )
              }
              
              
            } 
          })
        })
      } 
      
    })
  }
  )
}


# Run the application 
shinyApp(ui = ui, server = server)
