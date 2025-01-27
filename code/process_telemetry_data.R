# Mammals
library(sf)
library(terra)
library(ctmm)
library(sp)
library(spatialEco)

setwd("C:/evans/Croatia/mammals/data")
root = "C:/evans/Croatia/mammals"
  gpk  = file.path(root, "croatia_mammals.gpkg")
  
# croatia tz = "CEST"  
croatia.prj <- "EPSG:3765"
p = 0.4 # probability threshold for CDF

spp = c("bear", "wolf", "lynx")[1]

#*****************************************************
# Functions
UD.raster <- function(x, DF = c("CDF", "PMF"), mask.ud = TRUE, prj = NULL) {
  if(is.null(prj)) {
    proj <- attr(x,"info")$projection
  } else {
    proj <- prj
  }
  if(mask.ud) {
    ud.poly <- as.sf(UD, convex=FALSE)
	  sf::st_crs(ud.poly) <- sf::st_crs(proj) 
  }	
  DF = DF[1]
  dx <- x$dr[1]
  dy <- x$dr[2]
  xmn <- x$r$x[1] - dx / 2            
  xmx <- x$r$x[length(x$r$x)] + dx / 2
  ymn <- x$r$y[1] - dy / 2            
  ymx <- x$r$y[length(x$r$y)] + dy / 2										  
  e <- terra::ext(c(xmn, xmx, ymn, ymx))
  z <- x$r$z
  if(DF == "PMF") { 
    x <- x[["PDF"]] * prod(x$dr) 
  } else { 
    x <- x[[DF]] 
  }
  if(length(dim(x)) == 2) {
    x <- t(x[,dim(x)[2]:1])
  } else {
    x <- aperm(UD[,dim(UD)[2]:1,],c(2,1,3))
  }
  R <- terra::rast(x, crs = terra::crs(proj), extent = e)
    if(mask.ud) R <- mask(R, ud.poly)	  
  return(R)
}

#*****************************************************
# Bear individual autocorrelated kde fit
bear.opp <- st_read("bear_filtered_haphazard.shp") 
  bear.opp$source <- "opportunistic"
  bear.opp$CDF <- NA
    bear.opp <- bear.opp[,c("CDF", "source", "geometry")]
	
bear <- st_read(file.path(ddir, "bear", "Bear_GPS_All_HR.shp"))
iutm <- st_transform(bear, st_crs("EPSG:3767"))
  bear$timestamp <- as.POSIXct(as.character(paste(bear$UTC_DATE, bear$UTC_TIME)), 
             format="%Y-%m-%d %H:%M:%S", tz="UTC")
  names(bear)[which(names(bear) %in% c("CollarID", "LATITUDE", "LONGITUDE"))] <- c("individual.local.identifier", "location.long", "location.lat")    
bear <- bear[,c("individual.local.identifier", "timestamp", "location.long", "location.lat", "geometry")]

# remove duplicate time measurments and NA's
dx <- which(duplicated(bear$timestamp))
  if(length(dx) > 0) bear <- bear[-dx,]
na.idx <- unique(c(which(is.na(bear$individual.local.identifier))
                   which(is.na(bear$timestamp))))
  if(length(na.idx) > 0) bear <- bear[-na.idx,]

# table(bear$individual.local.identifier)

# Coerce to ctmm telemetry object
t.bear <- as.telemetry(st_drop_geometry(bear), timeformat="auto", timezone="UTC", 
                       projection=NULL, datum="WGS84")				   
  bear.kde <- list()
  bear.poly <- list()
    for(j in 1:length(t.bear)) {
      cat("Estimating KDE model for", names(t.bear)[j], "- obs", j, "of", length(t.bear), "\n") 
      i <- t.bear[[j]]
      projection(i) <- croatia.prj
        spp.sub <- bear[bear$individual.local.identifier == i@info$identity,]
          i@.Data[[5]] <- round(as.numeric(st_coordinates(spp.sub)[,1]) ,0) 
          i@.Data[[6]] <- round(as.numeric(st_coordinates(spp.sub)[,2]) ,0)
    	cat(nrow(spp.sub), "telemetry observsations", "\n")   
        mdl <- ctmm.guess(i, interactive=FALSE)
          mdl.ac <- ctmm.select(i, mdl)
            UD <- occurrence(i, CTMM = mdl.ac, variable = "utilization")
	      bear.poly[[j]] <- as.sf(UD, error=FALSE)
        rcdf <- UD.raster(UD, prj = croatia.prj) 
      bear.kde[[j]] <- raster.invert(rcdf)  
    }
    bear.pts <- do.call(rbind, lapply(bear.kde, \(i) { st_as_sf(as.points(i)) }))
      names(bear.pts)[1] <- "CDF"
      bear.pts$source <- "GPS_KDE"
	  
bear.pts <- bear.pts[bear.pts$CDF >= p,]
  plot(bear.pts["CDF"], pch=20)
    st_write(bear.pts, "croatia_mammals.gpkg", "bear") 

#*****************************************************
# Lynx group autocorrelated kde fit
lynx.opp <- st_read("lynx_filtered_haphazard.shp") 
  lynx.opp$source <- "opportunistic"
  lynx.opp$CDF <- NA
    lynx.opp <- lynx.opp[,c("CDF", "source", "geometry")]

lynx <- st_read(file.path(ddir, "lynx", "lynx_filtered_GPS.shp"))
  lynx$timestamp <- as.POSIXct(as.character(lynx$utc_date), format="%Y-%m-%d %H:%M:%S", tz="UTC")
    names(lynx)[which(names(lynx) %in% c("animal_id", "latitude", "longitude"))] <- c("location.long", "location.lat", "individual.local.identifier")    
lynx <- lynx[,c("individual.local.identifier", "timestamp", "location.long", "location.lat", "geometry")]

# remove duplicate time measurments
dx <- which(duplicated(lynx$timestamp))
  if(length(dx) > 0) lynx <- lynx[-dx,]

# table(lynx$individual.local.identifier)

t.lynx <- as.telemetry(st_drop_geometry(lynx), timeformat="auto", timezone="UTC", 
                       projection=NULL, datum="WGS84")
  lynx.kde <- list()
  lynx.poly <- list()
    for(j in 1:length(t.lynx)) {
      cat("Estimating KDE model for", names(t.lynx)[j], "- obs", j, "of", length(t.lynx), "\n") 
      i <- t.lynx[[j]]
      projection(i) <- croatia.prj
        spp.sub <- lynx[lynx$individual.local.identifier == i@info$identity,]
          i@.Data[[5]] <- round(as.numeric(st_coordinates(spp.sub)[,1]) ,0) 
          i@.Data[[6]] <- round(as.numeric(st_coordinates(spp.sub)[,2]) ,0)
    	cat(nrow(spp.sub), "telemetry observsations", "\n")   
        mdl <- ctmm.guess(i, interactive=FALSE)
          mdl.ac <- ctmm.select(i, mdl)
            UD <- occurrence(i, CTMM = mdl.ac, variable = "utilization")
        rcdf <- UD.raster(UD, prj = croatia.prj) 
          lynx.poly[[j]] <- as.sf(UD, error=FALSE)
      lynx.kde[[j]] <- raster.invert(rcdf)  
    }
    lynx.pts <- do.call(rbind, lapply(lynx.kde, \(i) { st_as_sf(as.points(i)) }))
      names(lynx.pts)[1] <- "CDF"
        lynx.pts$source <- "GPS_KDE"
    lynx.pts <- lynx.pts[lynx.pts$CDF >= p,]
      plot(lynx.pts["CDF"], pch=20)

lynx.pts <- rbind(lynx.pts, lynx.opp)
  st_write(lynx.pts, "croatia_mammals.gpkg", "lynx") 

#*****************************************************
# Wolf group autocorrelated kde fit
wolf.opp <- st_read("wolf_filtered_haphazard.shp") 
  wolf.opp$source <- "opportunistic"
  wolf.opp$CDF <- 0
    wolf.opp <-wolf.opp[,c("CDF", "source", "geometry")]

wolf <- st_read(file.path(ddir, "wolf", "wolf_filtered_GPS.shp"))
  wolf$timestamp <- as.POSIXct(as.character(wolf$datetime), format="%Y-%m-%d %H:%M:%S", tz="UTC")
    names(wolf)[which(names(wolf) %in% c("AnimalID", "LATITUDE", "LONGITUDE"))] <- c("location.lat", "location.long", "individual.local.identifier")    
wolf <- wolf[,c("individual.local.identifier", "timestamp", "location.long", "location.lat", "geometry")]

# remove duplicate time measurments
dx <- which(duplicated(wolf$timestamp))
  if(length(dx) > 0) wolf <- wolf[-dx,]

t.wolf <- as.telemetry(st_drop_geometry(wolf), timeformat="auto", timezone="UTC", 
                       projection=NULL, datum="WGS84")
  wolf.kde <- list()
  wolf.poly <- list()
    for(j in 1:length(t.wolf)) {
      cat("Estimating KDE model for", names(t.wolf)[j], "- obs", j, "of", length(t.wolf), "\n") 
      i <- t.wolf[[j]]
      projection(i) <- croatia.prj
        spp.sub <- wolf[wolf$individual.local.identifier == i@info$identity,]
          i@.Data[[5]] <- round(as.numeric(st_coordinates(spp.sub)[,1]) ,0) 
          i@.Data[[6]] <- round(as.numeric(st_coordinates(spp.sub)[,2]) ,0)
    	cat(nrow(spp.sub), "telemetry observsations", "\n")   
        mdl <- ctmm.guess(i, interactive=FALSE)
          mdl.ac <- ctmm.select(i, mdl)
            UD <- occurrence(i, CTMM = mdl.ac, variable = "utilization")
        rcdf <- UD.raster(UD, prj = croatia.prj) 
          wolf.poly[[j]] <- as.sf(UD, error=FALSE)
      wolf.kde[[j]] <- raster.invert(rcdf)  
    }
    wolf.pts <- do.call(rbind, lapply(wolf.kde, \(i) { st_as_sf(as.points(i)) }))
      names(wolf.pts)[1] <- "CDF"
        wolf.pts$source <- "GPS_KDE"
    wolf.pts <- wolf.pts[wolf.pts$CDF >= p,]
      plot(wolf.pts["CDF"], pch=20)

wolf.pts <- rbind(wolf.pts, wolf.opp)
  st_write(wolf.pts, "croatia_mammals.gpkg", "wolf") 


bio <- st_read("C:/evans/Croatia/data/Croatia.gpkg", "biogeographic_regions")
plot(wolf.pts["CDF"], pch=20)
  plot(st_geometry(bio), add=TRUE)

# birds
# Golden Eagle Aquila chrysaetos
# Eurasian griffon vulture Gyps fulvus
library(sf)
library(raster)
library(terra)
library(ctmm)
library(sp)
library(spatialEco)
library(SpatialKDE)

setwd("C:/evans/Croatia/birds")

# croatia tz = "CEST"  
croatia.prj <- "EPSG:3765"
p = 0.4 # probability threshold for CDF
bio <- st_read("C:/evans/Croatia/data/Croatia.gpkg", "biogeographic_regions")
mdl.spp = c("Aquila chrysaetos", "Gyps fulvus")
ref <- rast("C:/evans/Croatia/data/mask100m.tif")
kde.type = c("autocorrelated", "quartic")

#*****************************************************
# Functions
UD.raster <- function(x, DF = c("CDF", "PMF"), mask.ud = TRUE, prj = NULL) {
  if(is.null(prj)) {
    proj <- attr(x,"info")$projection
  } else {
    proj <- prj
  }
  if(mask.ud) {
    ud.poly <- as.sf(UD, convex=FALSE)
	  sf::st_crs(ud.poly) <- sf::st_crs(proj) 
  }	
  DF = DF[1]
  dx <- x$dr[1]
  dy <- x$dr[2]
  xmn <- x$r$x[1] - dx / 2            
  xmx <- x$r$x[length(x$r$x)] + dx / 2
  ymn <- x$r$y[1] - dy / 2            
  ymx <- x$r$y[length(x$r$y)] + dy / 2										  
  e <- terra::ext(c(xmn, xmx, ymn, ymx))
  z <- x$r$z
  if(DF == "PMF") { 
    x <- x[["PDF"]] * prod(x$dr) 
  } else { 
    x <- x[[DF]] 
  }
  if(length(dim(x)) == 2) {
    x <- t(x[,dim(x)[2]:1])
  } else {
    x <- aperm(UD[,dim(UD)[2]:1,],c(2,1,3))
  }
  R <- terra::rast(x, crs = terra::crs(proj), extent = e)
    if(mask.ud) R <- mask(R, ud.poly)	  
  return(R)
}

bw.Stoyan <- function(X, co = 0.15) {
    stopifnot(spatstat.geom::is.ppp(X))
    n <- spatstat.geom::npoints(X)
    W <- spatstat.geom::as.owin(X)
    a <- spatstat.geom::area.owin(W)
    stoyan <- co/sqrt(5 * n/a)
    return(stoyan)
}

#*****************************************************
# Golden Eagle "Aquila chrysaetos" individual 
# autocorrelated kde fit
eagle.opp <- st_read("C:/evans/Croatia/data/croatia_birds.gpkg", "Aquila chrysaetos")
  st_geometry(eagle.opp) <- "geometry" 
  eagle.opp$source <- "opportunistic"
  eagle.opp$CDF <- NA
    eagle.opp <- eagle.opp[,c("CDF", "source", "geometry")]

eagle.nests <- st_read(file.path(getwd(), mdl.spp[1], "data", "golden_eagle_nests.shp")) 
  st_geometry(eagle.nests) <- "geometry" 
  eagle.nests$source <- "nests"
  eagle.nests$CDF <- NA
    eagle.nests <- eagle.opp[,c("CDF", "source", "geometry")]	
	
# Read individual GPS files
eagle.gps <- lapply(c("Aq_chry_Surkan.shp", "Aqchr_Maleni.shp",
    "Aqchry_Dijana_aka_SuriSenj_052021_htrs.shp", "Aqchy_Oto_akaSurkan_htrs.shp"), 
    \(i)  { 
	  s <- st_read(file.path(getwd(), mdl.spp[1], "data", i)) 
	    if(!st_crs(s) == st_crs(croatia.prj)) { 
	      s <- st_transform(s, croatia.prj)
		}
     return(s)	
	})  	
gps <- eagle.gps[[1]] 
  gps$timestamp <- as.POSIXct(as.character(gps$GPSTime), dformat="%Y-%m-%d %H:%M:%S", tz="UTC")
    gps <- cbind(gps, st_coordinates(st_transform(gps, "EPSG:4326")))						  
	  names(gps)[which(names(gps) %in% c("GpsNmbr", "X", "Y"))] <- c("individual.local.identifier", "location.long", "location.lat")  	  
eagle.gps[[1]]  <- gps[,c("individual.local.identifier", "timestamp", "location.long", "location.lat", "geometry")]

gps <- eagle.gps[[2]] 
  gps$timestamp <- as.POSIXct(as.character(gps$GPSTime), dformat="%Y-%m-%d %H:%M:%S", tz="UTC")
    gps <- cbind(gps, st_coordinates(st_transform(gps, "EPSG:4326")))						  
	  names(gps)[which(names(gps) %in% c("GpsNmbr", "X", "Y"))] <- c("individual.local.identifier", "location.long", "location.lat")    
eagle.gps[[2]]  <- gps[,c("individual.local.identifier", "timestamp", "location.long", "location.lat", "geometry")]

gps <- eagle.gps[[3]] 
  gps$timestamp <- as.POSIXct(as.character(gps$UTC_dateti), dformat="%Y-%m-%d %H:%M:%S", tz="UTC")
    gps <- cbind(gps, st_coordinates(st_transform(gps, "EPSG:4326")))						  
	  names(gps)[which(names(gps) %in% c("device_id", "X", "Y"))] <- c("individual.local.identifier", "location.long", "location.lat")    
eagle.gps[[3]]  <- gps[,c("individual.local.identifier", "timestamp", "location.long", "location.lat", "geometry")]

gps <- eagle.gps[[4]] 
  gps$timestamp <- as.POSIXct(as.character(gps$UTC_dateti), dformat="%Y-%m-%d %H:%M:%S", tz="UTC")
    gps <- cbind(gps, st_coordinates(st_transform(gps, "EPSG:4326")))						  
	  names(gps)[which(names(gps) %in% c("device_id", "X", "Y"))] <- c("individual.local.identifier", "location.long", "location.lat")    
eagle.gps[[4]]  <- gps[,c("individual.local.identifier", "timestamp", "location.long", "location.lat", "geometry")]

eagle <- do.call(rbind, eagle.gps)
  #eagle <- eagle[which(lengths(st_intersects(eagle, bio)) > 0),]  
    unique(format(eagle$timestamp,"%Y"))

# remove duplicate time measurments and NA's
dx <- which(duplicated(eagle$timestamp))
  if(length(dx) > 0) eagle <- eagle[-dx,]
na.idx <- unique(c(which(is.na(eagle$individual.local.identifier)),
                   which(is.na(eagle$timestamp))))
  if(length(na.idx) > 0) eagle <- eagle[-na.idx,]

# table(eagle$individual.local.identifier)

if(kde.type == "autocorrelated") {
  # Coerce to ctmm telemetry object
  t.eagle <- as.telemetry(st_drop_geometry(eagle), timeformat="auto", timezone="UTC", 
                          projection=NULL, datum="WGS84")				   
  eagle.kde <- list()
  eagle.poly <- list()
    for(j in (1:length(t.eagle))) {
      cat("Estimating KDE model for", names(t.eagle)[j], "- obs", j, "of", length(t.eagle), "\n") 
      i <- t.eagle[[j]]
      projection(i) <- croatia.prj
        spp.sub <- eagle[eagle$individual.local.identifier == i@info$identity,]
          i@.Data[[5]] <- round(as.numeric(st_coordinates(spp.sub)[,1]) ,0) 
          i@.Data[[6]] <- round(as.numeric(st_coordinates(spp.sub)[,2]) ,0)
    	cat(nrow(spp.sub), "telemetry observations", "\n")   
        mdl <- ctmm.guess(i, interactive=FALSE)
          mdl.ac <- ctmm.select(i, mdl)
            UD <- occurrence(i, CTMM = mdl.ac, variable = "utilization")
	      eagle.poly[[j]] <- as.sf(UD, error=FALSE)
        rcdf <- UD.raster(UD, prj = croatia.prj) 
      eagle.kde[[j]] <- raster.invert(rcdf)  
    }
    eagle.pts <- do.call(rbind, lapply(eagle.kde, \(i) { st_as_sf(as.points(i)) }))
      names(eagle.pts)[1] <- "CDF"
      eagle.pts$source <- "GPS_KDE"	  
  eagle.pts <- eagle.pts[eagle.pts$CDF >= p,]
} else {
  win <- spatstat.geom::convexhull.xy(st_coordinates(eagle))
    x.ppp <- suppressWarnings(
      spatstat.geom::as.ppp(sf::st_coordinates(eagle)[,1:2], win))
  bw <- bw.Stoyan(x.ppp)
  khat <- rast(kde(eagle, band_width = bw, kernel = "triweight", 
               grid = raster(crop(ref, ext(eagle) )))

    as.points()
}

eagle.pts <- rbind(eagle.pts, eagle.opp)
  eagle.pts <- rbind(eagle.pts, eagle.nests)  
    st_write(eagle.pts, "croatia_birds.gpkg", mdl.spp[1]) 


st_write(eagle, "telemetry_bird_data.gpkg", "Aquila chrysaetos telemetry") 
st_write(eagle.opp, "telemetry_bird_data.gpkg", "Aquila chrysaetos opportunistic") 
st_write(eagle.nests, "telemetry_bird_data.gpkg", "Aquila chrysaetos nests") 


#*****************************************************
# Eurasian griffon vulture "Gyps fulvus" individual 
# autocorrelated kde fit
vulture.opp <- st_read("C:/evans/Croatia/data/croatia_birds.gpkg", "Gyps fulvus")
  st_geometry(vulture.opp) <- "geometry" 
  vulture.opp$source <- "opportunistic"
  vulture.opp$CDF <- NA
    vulture.opp <- vulture.opp[,c("CDF", "source", "geometry")]
vulture.feeding <- st_read(file.path(getwd(), mdl.spp[2], "data", "hraniliste_supovi_plan_Velebit.shp")) 
  st_geometry(vulture.feeding) <- "geometry" 
  vulture.feeding$source <- "foraging"
  vulture.feeding$CDF <- NA
    vulture.feeding <- vulture.feeding[,c("CDF", "source", "geometry")]	
	
# Read individual GPS files
vulture.gps <- read.csv(file.path(getwd(), mdl.spp[2], "data", "cleaned-data.csv"))
  vulture.gps <- st_as_sf(vulture.gps, coords = c("longitude", "latitude"), crs = "EPSG:3035", agr = "constant")
    vulture.gps <- st_transform(vulture.gps, croatia.prj) 
#vulture.gps <- vulture.gps[which(lengths(st_intersects(vulture.gps, bio)) > 0),]
  vulture.gps$timestamp <- as.POSIXct(as.character(vulture.gps$time), format="%m/%d/%Y %H:%M", tz="UTC")
    vulture.gps <- cbind(vulture.gps, st_coordinates(st_transform(vulture.gps, "EPSG:4326")))						  
	  names(vulture.gps)[which(names(vulture.gps) %in% c("name", "X", "Y"))] <- c("individual.local.identifier", "location.long", "location.lat")  	  
vulture <- vulture.gps[,c("individual.local.identifier", "timestamp", "location.long", "location.lat", "geometry")]

# remove duplicate time measurments and NA's
dx <- which(duplicated(vulture$timestamp))
  if(length(dx) > 0) vulture <- vulture[-dx,]
na.idx <- unique(c(which(is.na(vulture$individual.local.identifier)),
                   which(is.na(vulture$timestamp))))
  if(length(na.idx) > 0) vulture <- vulture[-na.idx,]

if(kde.type == "autocorrelated") {
  table(vulture.gps$individual.local.identifier)
  t.vulture <- as.telemetry(st_drop_geometry(vulture), timeformat="auto", timezone="UTC", 
                            projection=NULL, datum="WGS84")				   
  vulture.kde <- list()
  vulture.poly <- list()
    for(j in (1:length(t.vulture))) {
      cat("Estimating KDE model for", names(t.vulture)[j], "- obs", j, "of", length(t.vulture), "\n") 
	    flush.console()	 
        Sys.sleep(0.01)	
      i <- t.vulture[[j]]
      projection(i) <- croatia.prj
        spp.sub <- vulture[vulture$individual.local.identifier == i@info$identity,]
          i@.Data[[5]] <- round(as.numeric(st_coordinates(spp.sub)[,1]) ,0) 
          i@.Data[[6]] <- round(as.numeric(st_coordinates(spp.sub)[,2]) ,0)
    	cat(nrow(spp.sub), "telemetry observations", "\n")   
        mdl <- ctmm.guess(i, interactive=FALSE)
          mdl.ac <- ctmm.select(i, mdl)
            UD <- occurrence(i, CTMM = mdl.ac, variable = "utilization")
	      vulture.poly[[j]] <- as.sf(UD, error=FALSE)
        rcdf <- UD.raster(UD, prj = croatia.prj) 
      vulture.kde[[j]] <- raster.invert(rcdf)  
    }
    vulture.pts <- do.call(rbind, lapply(vulture.kde, \(i) { st_as_sf(as.points(i)) }))
      names(vulture.pts)[1] <- "CDF"
      vulture.pts$source <- "GPS_KDE"  
        vulture.pts <- vulture.pts[vulture.pts$CDF >= p,]
  } else {
    win <- spatstat.geom::convexhull.xy(st_coordinates(vulture))
      x.ppp <- suppressWarnings(
        spatstat.geom::as.ppp(sf::st_coordinates(vulture)[,1:2], win))
    bw <- bw.Stoyan(x.ppp)
    khat <- kde(vulture, band_width = bw, kernel = "triweight", 
                 grid = raster(crop(ref, ext(vulture) )))
    
  }  		
      vulture.pts <- rbind(vulture.pts, vulture.opp)
    vulture.pts <- rbind(vulture.pts, vulture.feeding)  
  st_write(vulture.pts, "croatia_birds.gpkg", mdl.spp[2]) 

# plot(vulture.pts["CDF"], pch=20)
#   plot(st_geometry(bio), add=TRUE)

st_write(vulture, "telemetry_bird_data.gpkg", "Gyps fulvus telemetry") 
st_write(vulture.opp, "telemetry_bird_data.gpkg", "Gyps fulvus opportunistic") 
st_write(vulture.feeding, "telemetry_bird_data.gpkg", "Gyps fulvus foraging") 
