invisible(lapply(c("sf", "geodata", "spatialEco"), 
          require, character.only=TRUE))

wgs84 <- "EPSG:4326"
etrs <- "EPSG:3035"

setwd("C:/evans/Croatia" )
  dat.dir <- file.path(getwd(), "data")

#*********************************************************
# Country boundary
#bdy <- st_as_sf(gadm("Croatia", level=0, path=tempdir(), 
#                version="latest", resolution=1))
#  geo.ext <- ext(bdy)
#    bdy <- st_transform(bdy, etrs)
# st_write(bdy, file.path(dat.dir, "bdy.shp"))

bdy <- st_read(file.path(dat.dir, "bdy.shp"))
  geo.ext <- ext(13.4895820620001, 19.4351768500001, 42.3854260000001, 46.550518036)

#*********************************************************
# 100m reference raster in ETRS projection
ref <- rast(ext(bdy), resolution = c(100,100), crs=crs(etrs)) 

#*********************************************************
# Create X,Y coordinates rasters

xcoord <- ref
  xcoord[] <- xyFromCell(ref, 1:ncell(ref))[,1]
    xcoord <- mask(xcoord, vect(bdy)) 
	  names(xcoord) <- "xcoord"
        writeRaster(xcoord, file.path(dat.dir, "xcoord.tif")) 
ycoord <- ref
  ycoord[] <- xyFromCell(ref, 1:ncell(ref))[,2]
    ycoord <- mask(ycoord, vect(bdy)) 
	  names(ycoord) <- "ycoord"
        writeRaster(ycoord, file.path(dat.dir, "ycoord.tif")) 
	  
#*********************************************************
# Landcover 100m fractional cover; trees, grassland, shrubs, wetland, cropland, built)
trees <- landcover("trees", path=tempdir())
  trees <- mask(project(crop(trees, geo.ext), ref, method="bilinear"), vect(bdy)) 
    names(trees) <- "trees"
      writeRaster(trees, file.path(dat.dir, "trees.tif")) 
grass <- landcover("grassland", path=tempdir())
  grass <- mask(project(crop(grass, geo.ext), ref, method="bilinear"), vect(bdy)) 
    names(grass) <- "grass"
      writeRaster(grass, file.path(dat.dir, "grass.tif"))   
shrub <- landcover("shrubs", path=tempdir())
  shrub <- mask(project(crop(shrub, geo.ext), ref, method="bilinear"), vect(bdy)) 
    names(shrub) <- "shrub"
      writeRaster(shrub, file.path(dat.dir, "shrub.tif")) 
wetland <- landcover("wetland", path=tempdir())
  wetland <- mask(project(crop(wetland, geo.ext), ref, method="bilinear"), vect(bdy)) 
    names(wetland) <- "wetland"
      writeRaster(wetland, file.path(dat.dir, "wetland.tif")) 
crops <- landcover("cropland", path=tempdir())
  crops <- mask(project(crop(crops, geo.ext), ref, method="bilinear"), vect(bdy)) 
    names(crops) <- "crop"
    writeRaster(crops, file.path(dat.dir, "crops.tif")) 

#*********************************************************
# Elevation (100m or 30 arc sec)

elev <- elevation_30s(country="HRV", path=tempdir() )
  elev <- mask(project(crop(elev, geo.ext), ref, method="bilinear"), vect(bdy))
    names(elev) <- "elev"
      writeRaster(elev, file.path(dat.dir, "elev.tif"), 
	              overwrite=TRUE)  

# Heat Load Index
heat.load <- hli(elev)
  names(heat.load) <- "hli"
    writeRaster(heat.load, file.path(dat.dir, "hli.tif"), 
	            overwrite=TRUE)  

# Topographic Roughness  
rough <- tri(elev, exact=FALSE)
  names(rough) <- "tri"
    writeRaster(rough, file.path(dat.dir, "tri.tif"), 
	            overwrite=TRUE)  

# Surface Relief Ratio 3x3 matrix
rr <- srr(elev, 3) 
  names(rr) <- "srr"
    writeRaster(rr, file.path(dat.dir, "srr3.tif"), 
	            overwrite=TRUE)  

# topographic position
tp <- tpi(elev, 3)
  names(tp) <- "tpi"
    writeRaster(tp, file.path(dat.dir, "tpi3.tif"), 
	            overwrite=TRUE)  

# Slope * COS(Aspect)
sa <- terra::terrain(elev, v=c("slope", "aspect"), unit="degrees")
  scosa <- terra::lapp(c(sa[[1]], sa[[2]]), fun = sa.trans)
    names(scosa) <- "scosa"
      writeRaster(scosa, file.path(dat.dir, "scosa.tif"), 
	              overwrite=TRUE)  

#*********************************************************
# WorldClim - temp min, max and precip
tmin <- worldclim_country("Croatia", var="tmin", path=tempdir())
  tmin <- median(tmin, na.rm=TRUE) 
    tmin <- mask(project(crop(tmin, geo.ext), ref, method="bilinear"), vect(bdy)) 
	  names(tmin) <- "tmin"
        writeRaster(tmin, file.path(dat.dir, "tmin.tif"), 
	                overwrite=TRUE)    
tmax <- worldclim_country("Croatia", var="tmax", path=tempdir())
  tmax <- median(tmax, na.rm=TRUE) 
   tmax <- mask(project(crop(tmax, geo.ext), ref, method="bilinear"), vect(bdy)) 
      names(tmax) <- "tmax"
        writeRaster(tmax, file.path(dat.dir, "tmax.tif"), 
	                overwrite=TRUE)     
tdiff <- tmin - tmax
  names(tdiff) <- "tdiff"
    writeRaster(tdiff, file.path(dat.dir, "tdiff.tif"), 
	            overwrite=TRUE)     
precip <- worldclim_country("Croatia", var="prec", path=tempdir())
  precip <- median(precip, na.rm=TRUE)
    precip <- mask(project(crop(precip, geo.ext), ref, method="bilinear"), vect(bdy)) 
	  names(precip) <- "precip"
        writeRaster(precip, file.path(dat.dir, "precip.tif"), 
	                overwrite=TRUE)  

#( r <- rast(list.files(dat.dir, "tif$", full.names=TRUE)) )
#  names(r)
  