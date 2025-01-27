library(terra)
library(sf)
library(spatialEco)

setwd("C:/evans/croatia/data")
ddir <- file.path(getwd(), "data100")

bdy <- st_read("Croatia.gpkg", "boundary")
eco <- st_cast(st_read("Croatia.gpkg", "biogeographic_regions"), "POLYGON")
  eco <- eco[eco$short_name == "mediterranean",]
lastovo <- st_read("Croatia.gpkg", "Dubrovnik-Neretva")
  lastovo <- lastovo[6,]
  
parms <- c("climate_pca.tif", "crops.tif", "development_multiscale.tif", "dist_dev.tif", 
           "dist_forest.tif", "dist_lake.tif", "dist_roads.tif", "dist_stream.tif", 
		   "forest_cover.tif", "forest_heights.tif", "forest_pct.tif", 
		   "grasslands_multiscale.tif", "hli.tif", "LAI_trend.tif", "scosa.tif", "srr.tif", 
            "tpi.tif", "tri.tif", "water.tif", "xcoord.tif", "ycoord.tif")
r <- rast(file.path(ddir, parms))	   

r.eco <- mask(crop(r, eco), eco)  
r.lastovo <- mask(crop(r, lastovo), lastovo)  

pdf("range_of_variability.pdf", height=8.5, width=11)
  for(i in 1:nlyr(r)) {
    n <- names(r)[i]
    y <- data.frame(ID = "mediterranean", y = na.omit(r.eco[[i]][])[,1])
      y$grp.mean <- mean(y$y)
    x <-  data.frame(ID = "lastovo", y = na.omit(r.lastovo[[i]][])[,1])
      x$grp.mean <- mean(x$y)
    dat <- rbind(y,x)
    p <- ggplot(dat, aes(x=y, fill=ID)) +
                geom_density(alpha=0.4, bw="sj") +
      geom_vline(data=dat, aes(xintercept = grp.mean, color = y),
                 linetype="dashed") +
        ggtitle(paste("Group PDF for", n)) + 
          xlab("PDF") +
            ylab(n)
    print(p)
  }
dev.off()
