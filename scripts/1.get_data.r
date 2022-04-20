
##############
## SET UP
##############
library("httr")
# library("getPass")
library("jsonlite")
library("tidyverse")
library("data.table")

# 'httr': to send requests and get response from SMARTER-backend API; 
# 'jsonlite': to parse JSON output, which is the default format of the API response
# 'getPass': (not strictly required) it will prompt for credentials in order to not store them in our code

## Generate a JWT token with R
# you need to generate a JWT token in order to get full access to smarter metadata


##############
## CONFIG FILE
##############
args = commandArgs(trailingOnly=TRUE)
if (length(args) >= 1) {
  
  #loading the parameters
  source(args[1])
  # source("~/config.R")
  
} else {
  #this is the default configuration, used for development and debug
  writeLines('Using default config')
  
  #this dataframe should be always present in config files, and declared
  #as follows
  config = NULL
  config = rbind(config, data.frame(
    base_folder = 'landscape_genomics_training',
    base_url = "https://webserver.ibba.cnr.it",
    usernm = ***REMOVED***,
    passwd = ***REMOVED***,
    species = "Goat", # Sheep or Goat
    force_overwrite = FALSE
  ))
}

writeLines(' - current values of parameters')
print(paste("base URL:", config$base_url))
print(paste("base folder is:", config$base_folder))
print(paste("user name:", config$usernm))
print(paste("password:", config$passwd))
print(paste("selected species:", config$species))

### sourcing the file with support R functions
fname = file.path(config$base_folder, "scripts","support_functions.r")
source(fname)

writeLines(" - getting the token for authentication")
token <- get_smarter_token(base_url = config$base_url, username = config$usernm, password = config$passwd)

writeLines(" - get list of breeds for the desired species")
breeds <- get_smarter_breeds(bsurl = config$base_url, token = token, query = list(species = config$species))

writeLines(" - get list of samples for the desired breed")
samples <- get_smarter_samples(base_url = config$base_url, token = token, species = config$species, query = list(breed_code = "ANK"))

writeLines(" - write out list of samples to filter the Plink binary file")
fname = file.path(config$base_folder,"samples_to_select.tsv")
select(samples, c(breed_code, smarter_id)) %>% fwrite(file = fname, sep="\t", col.names = FALSE)

