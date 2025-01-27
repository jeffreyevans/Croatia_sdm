# Create spatial data from parsed species flat-file
# NCBI entrez API key e99b3235858d25a0dbb4afbbcd08e5337308
library(sf)
library(spatialEco)

prj <- "EPSG:3035" # ESA ETRS89-extended / LAEA Europe projection https://epsg.io/3035
prj.croatia <- "EPSG:3765" # Croatian Terrestrial Reference System"

setwd("C:/evans/Croatia/species")
  d <- read.csv("all_birds.csv")
  # d <- read.csv("cleaned_table_ALL_bats.csv")
    na.idx <- unique(c(which(is.na(d$x)), which(is.na(d$y))))	
      if(length(na.idx) > 0) d <- d[-na.idx,]
	# Query common species names	  
	scn = lapply(unique(d$species), \(i) { c(i, unlist(taxize::sci2comm(sci=i))) })	  
	  scn.df <- as.data.frame(do.call(rbind, scn))
	    names(scn.df) <- c("species", "common_name")
  d <- dplyr::left_join(scn.df, d, by = "species")

spp.cts <- as.data.frame(table(d$species))
   spp.cts <- spp.cts[-which(spp.cts$Freq <= 20),]
d <- d[which(d$species %in% spp.cts$Var1),]

spp <- lapply(unique(d$species), \(i) {
    dspp <- d[d$species == i,]
      dspp <- st_as_sf(dspp, coords = c("x", "y"), crs = prj.croatia, agr = "constant")
	    dspp <- sf::st_difference(dspp)
		  dspp <- st_transform(dspp, prj)
	st_write(dspp, "croatia_bird_occurance.gpkg", i, append=FALSE)
	return(dspp)
})
names(spp) <- unique(d$species) 

   spp.nni <- lapply(spp, \(i) {
			  if(nrow(i) > 4) {
			    nn = nni(i)
			  } else {
                nn <- rep(NA,5)
				  names(nn) <- c("NNI", "z.score", "p", "expected.mean.distance", "observed.mean.distance")
              }
			  return(nn)
			})
    spp.nni <- lapply(spp.nni, \(i) { round(unlist(i), 5) } )		
  spp.nni <- as.data.frame(do.call(rbind, spp.nni))   
spp <- data.frame(scn.df, counts=unlist(lapply(spp, nrow)), spp.nni)
  spp <- spp[,-c(5,6)]
    names(spp)[5:6] <- c("expected_dist", "observed_dist")

write.csv(spp, "bird_summaries.csv")