library(sf)
library(readxl)

wgs84 <- st_crs("EPSG:4326")
htrs96 <- st_crs("EPSG:3765")
etrs <- st_crs("EPSG:3035")

root = "C:/evans/Croatia" 
setwd(root)
  spp.dir <- file.path(root, "spp")  # spp spreadsheet dir
  mdl.dir <- file.path(root, "SDM")  # root dir to create model dirs

f <- list.files(spp.dir, pattern = ".xlsx$", 
           full.names = TRUE)

spp_codes <- vector()
spp_names <- vector()
  for(i in f) {
    d <- as.data.frame(read_excel(i))
      if(length(table(d[,"Validno ime svojte"])) > 1) {
  	    sn <- table(d$species) 
  	    d$species <- names(sn)[which.max(sn)] 
      }
    spp <- sub(" ", "_", unique(d[,"Validno ime svojte"])[1]) 
    spp_code <- unique(tolower(gsub("^(.{3,4}).* (.{3,4}).*$", "\\1_\\2", 
                   d[,"Validno ime svojte"]))) 
    spp_names <- append(spp_names, spp)
    spp_codes <- append(spp_codes, spp_code)
    obs.type <- ifelse(as.numeric(d[,"Broj opaženog"]) == 0, "true abs", 
  	            ifelse(as.numeric(d[,"Broj opaženog"]) >= 1, "pres", NA))
    d <- data.frame(y = 1, obs_type = obs.type, species = d[,"Validno ime svojte"],
                    nobs = as.numeric(d[,"Broj opaženog"]),
  				  dttm = as.Date(paste0(d[,"Godina opažanja"], "-", 
  				  d[,"Mjesec opažanja"], "-", 
  				  d[,"Dan opažanja"]), "%Y-%m-%d"),
  				  xcoord = d[,"X koordinata"], 
  				  ycoord = d[,"Y koordinata"])
    d$y <- ifelse(d$nobs == 0, 0, d$y) 
    dsf <- st_transform(st_as_sf(d, coords = c("xcoord", "ycoord"), 
                        crs = htrs96), etrs)    
    dir.create(file.path(mdl.dir, spp), showWarnings = FALSE)			 
    st_write(dsf, file.path(mdl.dir, spp, paste0(spp,".shp")), 
	         quiet = TRUE, delete_layer = TRUE)	
	cat("\n")
    cat("spp point data written to:", file.path(mdl.dir, spp, paste0(spp,".shp")), "\n")
    cat("\n")	
  }  
