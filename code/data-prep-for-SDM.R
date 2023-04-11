# 
# SMART ciljevi - data analysis
# Data prep for Random Forest modelling
# Biom 2023
# 

# uncomment to install required packages (also add if any are missing):
# install.packages(c("tidyverse", "readxl", "sf"))

# load packages:
library(tidyverse)
library(readxl)
library(sf)

wgs84 <- st_crs("EPSG:4326")
htrs96 <- st_crs("EPSG:3765")
etrs <- st_crs("EPSG:3035")

vect_of_input_tables <- list.files(path = "data-input/observations/example_data", 
                                   pattern = ".xlsx", 
                                   full.names = TRUE)

get_select_transform <- function(filename) {
  output <- read_excel(filename) %>% 
    st_as_sf(coords = c("X koordinata", "Y koordinata"), crs = htrs96) %>% 
    st_transform(crs = etrs)
  
  # drop unnecessary columns
  output <- output %>% transmute(species = `Validno ime svojte`,
                           n_obs = as.numeric(`Broj opaženog`),
                           dttm = ymd(paste0(`Godina opažanja`, "-", `Mjesec opažanja`, "-", `Dan opažanja`)))
  return(output)
}

tab <- lapply(vect_of_input_tables, FUN = "get_select_transform") %>% 
  bind_rows()

# check for na values in species and observations of multiple birds
tab %>% filter(is.na(species))
tab %>% filter(n_obs > 1)

# check for typos and missing values in species column
obs_species <- tab$species %>% unique
obs_species

# for each transect, need list of birds that would have been recorded
# if they were present in the area -- to generate absences where appropriate

for (i in obs_species) {
  # get species code in format Aaaa_bbbb
  speccode <- gsub("^(.{3,4}).* (.{3,4}).*$", "\\1_\\2", i)
  outputdir <- paste0("data-output/", speccode)
  try(dir.create(outputdir), silent = TRUE)
  
  tab_i <- filter(tab, species == i)
  
  # write out presence data
  tab_i %>% filter(n_obs > 0) %>% 
    write_sf(dsn = paste0(outputdir,"\\", speccode, "_pres_etrs.shp"))

  # check for absences and write them out if they exist
  if(any(tab_i$n_obs == 0)) {
    tab_i %>% filter(n_obs == 0) %>% 
      write_sf(dsn = paste0(outputdir,"\\", speccode, "_abs_etrs.shp"))
  }
}

# check for na values in species and observations of multiple birds
tab %>% filter(is.na(species))
tab %>% filter(n_obs > 1)

# check for typos and missing values in species column
obs_species <- tab$species %>% unique
obs_species

# for each transect, need list of birds that would have been recorded
# if they were present in the area -- to generate absences where appropriate

try(dir.create("data-input/prepared-shapefiles"))

for (i in obs_species) {
  # get species code in format Aaaa_bbbb
  speccode <- gsub("^(.{3,4}).* (.{3,4}).*$", "\\1_\\2", i)
  outputdir <- paste0("data-input/prepared-shapefiles/", speccode)
  try(dir.create(outputdir), silent = TRUE)
  
  tab_i <- filter(tab, species == i)
  
  # write out presence data
  tab_i %>% filter(n_obs > 0) %>% 
    write_sf(dsn = paste0(outputdir,"/", speccode, "_pres_etrs.shp"))
  
  # check for absences and write them out if they exist
  if(any(tab_i$n_obs == 0)) {
    tab_i %>% filter(n_obs == 0) %>% 
      write_sf(dsn = paste0(outputdir,"/", speccode, "_abs_etrs.shp"))
  }
}



