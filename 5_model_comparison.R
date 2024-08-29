############################################################################################
################  collate model stats on dev exp, null dev, perc explained   ############### 
###############   and identify the best model to use in Reef Adapt           ############### 
############################################################################################
# Written by G V Wood and K.J. Griffin Feb 2022
# Email george.wood@flinders.edu.au for additional support

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("0_Dependencies.R")


# List of directories and their corresponding model names
directory_mapping <- list(
  "output/gdmmodel/connectivity_unscaled/" = "conn_unsc",
  "output/gdmmodel/connectivity_scaled/" = "conn_sc",
  "output/gdmmodel/connectivity_scaled_euc/" = "conn_sc_euc",
  "output/gdmmodel/connectivity_unscaled_euc/" = "conn_unsc_euc"
)

# Initialize an empty data frame to store results
gdm_results <- data.frame()

# Loop through each directory
for(dir in names(directory_mapping)) {
  
  # Get the model name corresponding to the directory
  model_name <- directory_mapping[[dir]]
  
  # Get the list of files in the current directory
  rds_combo <- list.files(path = dir, pattern = "modTest_", full.names = FALSE) 
  rds_combo_fullnames <- list.files(path = dir, pattern = "modTest_", full.names = TRUE) 
  
  # Loop through each file in the current directory
  for(i in 1:length(rds_combo)){
    skip_to_next <- FALSE
    
    tryCatch({
     
      # Extract refid
      get_refid <- function(rds_filename) {
        # Try to match with underscores first
        refid <- gsub(".*(ref_[0-9]+(?:_[a-z0-9]+)*)\\.rds", "\\1", rds_filename)
        
        # If no underscore match, try matching with numbers and optional trailing lowercase letters
        if (refid == rds_filename) {
          refid <- gsub(".*(ref[0-9]+(?:_[a-z]+)?)\\.rds", "\\1", rds_filename)
        }
        
        return(refid)
      }

      refid <- get_refid(rds_combo[i])
      
      # Find the corresponding bootstrap results file
      bootstrap_results_rds <- list.files(path = dir, pattern = paste0("modTest_", refid, ".rds"), full.names = TRUE) 
      
      # Read the bootstrap results
      bootstrap_results <- readRDS(bootstrap_results_rds)
      
      # Create a data frame with the results
      gdm_bootstrap_results_df <- data.frame(
        "refid" = refid,
        "model" = model_name,
        "gdmdev" = bootstrap_results[[1]][1,],
        "nulldev" = (bootstrap_results[[1]][1,])/(1-(bootstrap_results[[1]][2,]/100)),
        "explained" = bootstrap_results[[1]][2,],
        "pval" = bootstrap_results[[1]][3,],
        "fitted_perms" = bootstrap_results[[1]][4,],
        "geo_imp" = bootstrap_results[[2]][1,],
        "env_imp" = sum(bootstrap_results[[2]][-1,])
      )
      
      # Append the results to the main data frame
      gdm_results <- rbind(gdm_results, gdm_bootstrap_results_df)
      
    }, error = function(e){skip_to_next <- TRUE})
    
    # Skip to the next iteration if an error occurred
    if(skip_to_next){next}
  }
}

gdm_results
unique(gdm_results$refid)
############################################################################################
################               Curate which models should go forward         ############### 
############################################################################################
# Create a dataframe to store the user's choices
chosen_models <- data.frame(refid = character(), model = character(), stringsAsFactors = FALSE)

# Get unique refids
unique_refids <- unique(gdm_results$refid)

# Define the model choices
model_choices <- c("conn_unsc", "conn_sc", "conn_sc_euc", "conn_unsc_euc")

for (refid in unique_refids) {
  # Filter the results for the current refid
  refid_results <- gdm_results[gdm_results$refid == refid, ]
  
  # Check if the top explained value is less than 20
  max_explained <- max(refid_results$explained)
  if (max_explained < 20) {
    print(paste("Top explained value for refid", refid, "is less than 20. Automatically excluding this refid."))
    next
  }
  
  # Display the results for the current refid
  print(paste("Results for refid:", refid))
  print(refid_results)
  
  # Display model choices
  print("Model choices:")
  print("0 : Exclude this refid")
  for (i in 1:length(model_choices)) {
    print(paste(i, ":", model_choices[i]))
  }
  
  # Prompt the user to choose the best model by number
  choice <- as.integer(readline(prompt = "Enter the number corresponding to the model you wish to go ahead with (0 to exclude): "))
  
  # Validate choice
  while (is.na(choice) || choice < 0 || choice > 4) {
    choice <- as.integer(readline(prompt = "Invalid choice. Please enter a number between 0 and 4: "))
  }
  
  # Store the user's choice in the chosen_models dataframe
  if (choice != 0) {
    chosen_model <- model_choices[choice]
    chosen_models <- rbind(chosen_models, data.frame(refid = refid, model = chosen_model, stringsAsFactors = FALSE))
  }
}

###### PAUSE HERE AND CHOOSE FINAL MODELS ######

# Display the final chosen models
print(chosen_models)

# Check that these models are oK to go ahead with:
# Subset the gdm_results dataframe to only include the chosen models
final_gdm_results <- merge(gdm_results, chosen_models, by = c("refid", "model"))

# Save the final chosen models to a file 
dir.create("output/gdm_model_output/useful_models/")
write.csv(final_gdm_results, "output/gdm_model_output/useful_models/chosen_models.csv", row.names = FALSE)




