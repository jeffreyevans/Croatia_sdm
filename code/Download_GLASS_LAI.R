# Download GLASS LAI (250m) or Fractional Vegetation Cover (500m)
# in 8 day composites
#
#   http://www.glass.umd.edu/introduction.html
#
# Jeffrey S. Evans, Ph.D.,
# Senior Landscape Ecologist & Biometrician 
# The Nature Conservancy | Global Protect, Science 
# jeffrey_evans@tnc.org 
#
#**************************************************
# Add libraries
suppressMessages(
  lapply(c("sf", "spatialEco", "terra", "rvest", "stringr"), 
         require, character.only = TRUE))

#**************************************************
# Define product and working directory
product.idx = 1 # 1 = LAI (250m), 2 = Fractional Vegetation Cover (500m)
url = c("http://www.glass.umd.edu/LAI/MODIS/250m",
        "http://www.glass.umd.edu/FVC/MODIS/500m/")[product.idx]
out.res = c(250, 500)[product.idx]

# Croatia ETRS89-extended / LAEA Europe projection https://epsg.io/3035
prj <- "EPSG:3035"

#***********************************************************
# Read gpkg study area boundary
bdy <- st_read(file.path(ddir, "Croatia.gpkg"), "boundary")
ref <- rast(ext(bdy), resolution = 30, crs=crs(prj))
 	
dat.dir = "C:/evans/Croatia/data/timeseries" #*** User defined 
  setwd(dat.dir)

out.file <- "LAI_250m_raw.tif" #*** User defined

#**************************************************
# Add vector polygon that defines study area and create 
# reference raster
#
# read vector polygon data representing download area. 
#   note; you could create a polygon bounding box from 
#         data and even buffer it to define padded area
#     eg., bdy <- st_as_sf(as.polygons(ext(bdy)), crs = st_crs(bdy) )
#            bdy <- st_buffer(bdy, 2500) 
bdy <- st_read("C:/evans/Croatia/data/Croatia.gpkg", "boundary") #*** User defined
  ref <- rast(vect(bdy), resolution = out.res, crs = crs(bdy))
    ref <- rasterize(bdy, ref)

# read MODIS tile grid
#*** path is user defined
tiles <- st_read(file.path("C:/evans/remote_sensing/tile_grids/MODIS",
                 "modis_sinusoidal_grid_world.shp"))
    bdy.prj <- st_transform(bdy, st_crs(tiles))

# Create tile index from intersecting study area intersection
ids <- st_drop_geometry(st_intersection(tiles, st_transform(bdy, st_crs(tiles))))[,c("h","v")]
  h <- unique(ids[,"h"])
    h <- str_pad(h, 2, pad = "0")
  v <- unique(ids[,"v"])
    v <- str_pad(v, 2, pad = "0")
modis.tiles <- do.call(paste, c(expand.grid("h", h, "v", v), sep = ""))
#modis.tiles <- c("h10v03", "h09v03")

# Define years (available 2000-2021)
year.range <- 2000:2021 #*** User defined
  year.range <- split(year.range, ceiling(seq_along(year.range)/6))[[4]]

#**************************************************
# Create year directory index and url's
sess <- session(url)
  years <- as.data.frame(html_table(sess)[[1]][,2])[,1]
    years <- years[3:24]
      years <- grep(paste(year.range, collapse="|"), years, value = TRUE) 
d <- lapply(years, \(i) {
 s <- session(paste0(url, "/", i))
  tiles <- as.data.frame(html_table(s)[[1]][,2])[,1]
    tiles <- tiles[-c(1,2,length(tiles))]
    paste0(url, "/", i, tiles)
})
names(d) <- years

#**************************************************
# Parse Julian day directories and create URL's
# for specified tiles 
f <- lapply(d, \(i) {
   sapply(i , FUN = \(j) {
     hf <- as.data.frame(html_table(session(j))[[1]][,2])[,1]
       hf <- grep("hdf", hf, value = TRUE)
         hdf.files <- hf[-grep("xml", hf)]
        if(exists("modis.tiles")) {
          hdf.files <- grep(paste(modis.tiles, collapse="|"), hdf.files, value=TRUE)
        }         
      return(paste0(j, hdf.files))
    }, simplify = FALSE, USE.NAMES = TRUE)
})
names(f) <- years

#**************************************************range
# Connect to rasters, mosaic, project and mask to 
# defined study area. Please note that the Julian
# day is converted to a data and is defined as the
# rasters name, which can be parsed later
for(i in 1:length(f)) {
  cat("Processing", i, "in", length(f), "for", names(f)[i], "\n") 
  y <- f[[i]]
  n <- unlist(lapply(strsplit(names(y), "/"), \(p) paste0(p[8], "_", p[7])))
    jdate <- as.Date(unlist(lapply(strsplit(n, "_"), \(y) {
      as.Date(paste(y[2], y[1]), format="%Y %j")
    })))
  for(p in 1:length(y)) {
    if(!file.exists(file.path(dat.dir, paste0("lai_", n[p], ".tif")))){  
      cat("Mosaic and projecting", p, "in", length(y), "for", n[p], "\n")
	  #r <- sprc(lapply(y[[p]], \(i) crop(rast(i), ext(bdy.prj))))
      #  r <- mosaic(r)
      #    r <- project(r, ref)
	  r <- sprc(lapply(y[[p]], \(i) project(rast(i), ref)))
        r <- mosaic(r)
          r <- mask(ifel(r > 9, NA, r), vect(bdy))
            names(r) <- jdate[p]
          writeRaster(r, file.path(dat.dir, paste0("lai_", n[p], ".tif")),
                      overwrite = TRUE, gdal=c("COMPRESS=LZW"), 
                      datatype="FLT4S")
		suppressWarnings({			  
	      try(lapply(file.path("C:/temp", basename(y[[p]])),  file.remove))
          terra::tmpFiles(current=TRUE, orphan=TRUE, old=TRUE, remove=TRUE)
		})	
      gc()		   
    } else {
      cat(paste0("lai_", n[p], ".tif"), "already exists", "\n")     
    }
  } 
}

#**************************************************
# Create raster timeseries from individual files
#
# LAI tif files
f <- intersect(list.files(dat.dir, pattern = "lai", full.names=TRUE), 
               list.files(dat.dir, pattern = "tif$", full.names=TRUE))

# Create dates and sorting index from file names
dates <- as.Date(unlist(lapply(strsplit(rm.ext(basename(f)), "_"), \(i) {
  as.Date(paste(i[3], i[2]), format="%Y %j")
 })))
  date.idx <- sort.int(dates, index.return=TRUE)$ix

# read rasters in order of timeseries
r <- rast(f[date.idx])
  names(r) <- dates[date.idx]

writeRaster(r, file.path(dirname(getwd()), out.file), 
            overwrite = TRUE, 
            gdal=c("COMPRESS=LZW"), 
            datatype="FLT4S")
