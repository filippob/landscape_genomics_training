
##############
## SET UP
##############

library(smarterapi)
library(dplyr)
library(data.table)

################
## CONFIGURATION
################
config <- list(
  base_folder = 'landscape_genomics_training',
  species = "Goat" # Sheep or Goat
)

writeLines(' - current values of parameters')
print(paste("base folder is:", config$base_folder))
print(paste("selected species:", config$species))

# token is managed through smarterapi

writeLines(" - get list of breeds for the desired species")
breeds <- smarterapi::get_smarter_breeds(query = list(species = config$species))

writeLines(" - get list of samples for the desired breed")
samples <- smarterapi::get_smarter_samples(species = config$species, query = list(breed_code = "ANK"))

writeLines(" - write out list of samples to filter the Plink binary file")
fname = file.path("selected_samples.tsv")
dplyr::select(samples, c(breed_code, smarter_id)) %>% data.table::fwrite(file = fname, sep = "\t", col.names = FALSE)

