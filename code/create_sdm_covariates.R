invisible(lapply(c("sf", "geodata", "spatialEco"), 
          require, character.only=TRUE))

wgs84 <- "EPSG:4326"
etrs <- "EPSG:3035"

setwd("C:/evans/Croatia")
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

#***********************************************************
# Create STAC URL for Copernius DEM digital surface model COG
ref <- rast(ext(bdy), resolution = 30, crs=crs(prj))
 
if(!file.exists(file.path(ddir, "elev.tif"))) {
  olm <- read_stac("http://s3.eu-central-1.wasabisys.com/stac/openlandmap/catalog.json")
    olm$links <- links(olm, rel == "child")
      links(olm, grepl("DEM", title))
  glc_link <- links(olm, grepl("DEM", title))[[1]]
    glc <- link_open(glc_link)
  glc_items <- read_items(glc, progress = FALSE)
    items_assets(glc_items)
  urls <- assets_url(glc_items, asset_names = "dsm_glo30_m_30m_s", append_gdalvsi = TRUE)
  
  # crop and reproject DEM's from COG, write raster
  e <- ext(st_transform(bdy, st_crs(4326)))
  elev <- rast(urls) 
    elev <- crop(elev, e)
      elev <- mask(project(elev, ref, method = "bilinear"), bdy)
        names(elev) <- "elevation" 
  writeRaster(elev, file.path(ddir, "elev.tif"),
              overwrite = TRUE, gdal=c("COMPRESS=LZW"), 
              datatype="FLT4S")
} else {
  elev <- rast(file.path(ddir, "elev.tif"))
}

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

if(!file.exist(file.path(dat.dir, "elev_100m.tif")))
  elev <- resample(rast(file.path(dat.dir, "elev.tif")), ref)
    names(elev) <- "elev"
    writeRaster(elev, file.path(dat.dir, "elev_100m.tif"), 
	            overwrite=TRUE)  
}

# Heat Load Index
sa <- terra::terrain(elev, v=c("slope", "aspect"), unit="degrees")
hli2 <- hli.new(sa[["aspect"]], sa[["slope"]], equation = 2)
  names(hli2) <- "hli"
    writeRaster(hli2, file.path(ddir, "hli.tif"), 
	            overwrite = TRUE, gdal=c("COMPRESS=LZW"), 
                datatype="FLT4S")

# Topographic Roughness  
rough <- c(tri(elev, exact=FALSE),
           tri(elev, s=11, exact=FALSE),
           tri(elev, s=27, exact=FALSE))
  rough <- round(rough, 4)	   
  names(rough) <- c("tri3", "tri11", "tri27")
    writeRaster(rough, file.path(ddir, "tri.tif"), 
	            overwrite = TRUE, gdal=c("COMPRESS=LZW"), 
                datatype="FLT4S")

# Surface Relief Ratio 3x3 matrix
rr <- c(srr(elev, 3),
        srr(elev, 11),
        srr(elev, 27))
  rr <- round(rr, 4)	   
  names(rr) <- c("srr3", "srr11", "srr27")
    writeRaster(rr, file.path(ddir, "srr.tif"), 
	            overwrite = TRUE, gdal=c("COMPRESS=LZW"), 
                datatype="FLT4S")

# topographic position
tp <- c(tpi(elev, 3),
        tpi(elev, 11),
        tpi(elev, 27))
  tp <- round(tp, 4)	   		
  names(tp) <- c("tpi3", "tpi11", "tpi27")
    writeRaster(tp, file.path(ddir, "tpi.tif"), 
	            overwrite = TRUE, gdal=c("COMPRESS=LZW"), 
                datatype="FLT4S")

# Slope * COS(Aspect)
if(!exists("sa")) {
  sa <- terra::terrain(elev, v=c("slope", "aspect"), unit="degrees")
}
scosa <- terra::lapp(c(sa[[1]], sa[[2]]), fun = sa.trans)
  scosa <- round(scosa, 4)
  names(scosa) <- "scosa"
    writeRaster(scosa, file.path(ddir, "scosa.tif"), 
                overwrite = TRUE, gdal=c("COMPRESS=LZW"), 
                datatype="FLT4S")
				  

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
  