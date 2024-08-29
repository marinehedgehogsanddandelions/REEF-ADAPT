############################################################################################
############################################################################################

# This script loads the libraries and functions needed to build and analyse GDM models for Wood et al (unpub)

# Written by G V Wood and K.J. Griffin Feb 2022
# Email george.wood@flinders.edu.au for additional support

############################################################################################
###########################################################################################

# Package names
packages <- c("sdmpredictors","sf","gecko", "fasterize", "JuliaConnectoR","filesstrings", "vegan","terra","geodata","leaflet", "gdm", "raster", "sp", "formattable", "purrr" , "dplyr", "tidyverse", "data.table", "plotly", "ggplot2", "spaa", "stringr","tidyr","reshape2","graph4lg","caret", "dplyr", "ggplot2","ggplotlyExtra","plotly", "plyr", "gridExtra","stringr", "effectsize","RStoolbox", "htmltools", "leafem", "plainview","stars","mapview","leafsync","leaflet.extras2","leafpop","slideview","htmlwidgets", "ggpubr")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))

