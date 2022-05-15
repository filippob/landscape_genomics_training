
##############
## SET UP
##############
require("here")
require("dplyr")
require("data.table")
require("smarterapi")

workdir <- here("scripts")
message("Moving to the working directory: ", workdir)
setwd(workdir)

# token is managed through smarterapi
samples <- get_smarter_samples(
  species = "Goat",
  query = list(
    country="Italy",
    type="background",
    breed="Orobica",
    breed="Aspromontana",
    breed="Bionda dell'Adamello",
    breed="Argentata"
  )
)

writeLines(" - write out list of samples to filter the Plink binary file")
fname = file.path(here("selected_samples.tsv"))
dplyr::select(samples, c(breed_code, smarter_id)) %>% data.table::fwrite(file = fname, sep = "\t", col.names = FALSE)

