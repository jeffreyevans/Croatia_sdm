library(future)
library(future.apply)

invisible(lapply(c("sf", "spatialEco", "terra", "dplyr", 
          "parallel", "snow", "doSNOW"), 
           require, character.only = TRUE))

#***************************************
#***************************************
# set up environment
metrics <- "lai"
res <- c("250m", "300m", "1000m")[1]
country <- "croatia"

root = file.path("Z:/croatia", "data")
setwd(root)
  dat.dir = getwd()
  
autocorrelation = c(TRUE, FALSE)[1]
forest.mask = c(TRUE, FALSE)[2]
np = availableCores()-1

#***************************************
#***************************************
# functions
mk.tau <- function(x, min.obs = 10, offset = 5, fdate = c(2014, 1), freq = 36, 
            type = c("arima", "lm", "yuepilon", "zhang", "hamed", "none")) {
  knames <- c("tau", "trend", "slope", "intercept", "p.value", "LCL", "UCL")
  type = type[1]
  if(type %in% c("zhang", "yuepilon")) { 
    n = 7 
  } else if(type == "hamed") {
    n = 3   
  } else if(type %in% c("arima", "lm", "none")) {
    n = 6
  }
  if(length(x[!is.na(x)]) < min.obs) {
    tau <- rep(NA,n) 
    } else {
      if(type == "Zhang") {
  	    kcol <- c("tau", "trend", "trendp", "sig", "intercept", "lbound", "ubound")
  	      tau <- round(zyp::zyp.zhang(x[!is.na(x)]),5)[kcol]
  		    names(tau) <- knames 
      } else if(type == "hamed") {
        kcol <- c("Tau", "Sen's slope", "new P-value") 
          tau <- round(modifiedmk::mmkh(x[!is.na(x)],ci=0.95),5)[kcol] 
  		    names(tau) <- knames[c(1,3,5)] 
      }  else if(type == "yuepilon") {
  	    kcol <- c("tau", "trend", "trendp", "sig", "intercept", "lbound", "ubound")
  	      tau <- round(zyp::zyp.yuepilon(x[!is.na(x)]),5)[kcol]
  		    names(tau) <- knames
      }  else if(type == "arima") {
  	    kcol <- c("tau", "slope", "intercept", "p-value", "limits.LCL", "limits.UCL") 
        am <- tryCatch({
            x - as.numeric(resid(arima(ts(x, start=fdate, frequency = freq), order=c(1,0,0))) )
          }, error = function(e) {
            x - as.numeric(resid(forecast::auto.arima(ts(x, start=fdate, frequency = freq))))
        }) 
        tau <- round(spatialEco::kendall(am),5)[kcol] 
          names(tau) <- knames[-2]	 
      }  else if(type == "none") {
  	    kcol <- c("tau", "slope", "intercept", "p-value", "limits.LCL", "limits.UCL") 
          tau <- round(spatialEco::kendall(x[!is.na(x)]),5)[kcol] 
            names(tau) <- knames[-2]
      }  else if(type == "lm") {
        kcol <- c("tau", "slope", "intercept", "p-value", "limits.LCL", "limits.UCL")
          p <- predict(lm(x[-c(1:offset)]~x[-c((length(x)-(offset-1)):length(x))]))		
            tau <- round(spatialEco::kendall(p),5)[kcol]
			  names(tau) <- knames[-2]	 
	  }
  }
  return( tau )
}

#***************************************
#***************************************
# read raster time-series and parse dates
metric <- rast(file.path(root, paste0(toupper(metrics), "_", res, "_detrend", ".tif")))
  if(res == "250m") { 
   dates <- as.Date(names(metric))
  } else {
    d <- unlist(lapply(strsplit(names(metric), "_"), function(x) x[2]))
    dates <- as.Date(unlist(lapply(d, function(j) {
        paste( substr(j,1,4), substr(j,5,6), substr(j,7,8), sep="-" )
      })))
  }

if(forest.mask) {
  f <- rast(file.path(dat.dir, "forest_300m.tif"))
  metric <- mask(metric, f)
}

#***************************************
#***************************************
# Calculate Mann-Kendall TAU and Sen Slope
# x <- as.numeric(metric[475838]) # test
# mdls <- c("yuepilon", "zhang", "hamed", "arima", "lm", "none")
#   kr <- lapply(mdls, \(i) mk.tau(x, type=i) ) 
#     names(kr) <- mdls

#***********************
# without autocorrelation correction 
if(autocorrelation == FALSE) {
  metric.tau <- raster.kendall(metric)
} else {			  
#***********************
# Autocorrelation corrected with
# decorrelated timeseires using ARIMA
tau.dat <- as.data.frame(metric, cells=TRUE, na.rm=NA)
  cell.ids <- tau.dat[,1] 
    tau.dat <- tau.dat[,-1]
  
# Future multiprocessing
  plan("future::multisession", workers = np, gc = TRUE)
    system.time({
      tau <- future_apply(tau.dat, MARGIN=1, mk.tau, type = "arima")
	  #tau <- apply(tau.dat, MARGIN=1, mk.tau, type = "arima")
  })
  plan(sequential)

#cl <- makeCluster(parallel::detectCores()-1, type='SOCK')
#  clusterEvalQ(cl, c(library(tsfeatures)))
#    clusterExport(cl, c("process.ts"), 
#                  envir=environment())
#  registerDoSNOW(cl)
#    tse <- parApply(cl, tdat, MARGIN=1, FUN = function(o) { process.ts(na.omit(o)) }) 
#  stopCluster(cl)
#registerDoSEQ()

tau <- as.data.frame(t(tau))
  # set insignificant tau's and slope's to 0
  tau$tau <- ifelse(tau$p.value > 0.05, 0, tau$tau) 
  tau$slope <- ifelse(tau$p.value > 0.05, 0, tau$slope) 
metric.tau <- metric[[1]]
  metric.tau[] <- rep(NA, ncell(metric.tau))
    metric.tau <- rep(metric.tau, ncol(tau))
      names(metric.tau) <- names(tau)
metric.tau[cell.ids] <- tau

terra::writeRaster(metric.tau, file.path(mdl.dir, 
                   paste0(toupper(metrics), "_tau", ".tif")),
                   overwrite = TRUE)
}
				
			
