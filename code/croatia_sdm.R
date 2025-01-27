# Circaetus gallicus
invisible(lapply(c("sf", "spatialEco", "terra", "ranger", 
          "ggplot2", "randomForest", "rfUtilities", "pdp"), 
		  require, character.only=TRUE))

set.seed(42)                            # Random seed, for reproducibility 

group = c("birds", "mammals", "bats")[1]
spp.dat = paste0("croatia_", group, ".gpkg") 

root = "Z:/croatia"                     # Main root directory 
setwd(root)                             # Set working directory to main root directory
  dat.dir <- file.path(root, "data")    # directory with raster covariates
  mdl.dir <- file.path(root, group)     # root directory of models

# Croatian Terrestrial Reference System"
prj <- "EPSG:3765"

# number of pseudo-absence samples, if nsample = NULL  
# then the defined sample fraction (sfract) will be used, 
# based on number 0f non-NA raster, to draw sample size.
# At 100m are 237,067 no-NA cell values
nsample = 1000  # or NULL
sfract = 0.01   # eg., 1% sample of raster cells (n=2371)                       
nboot = 1001                                    # number of bootstrap replicates
bw.method = c("Diggle", "Scott")[1]             # Bandwidth plug-in (1st or 2nd order)
min.ct = 20                                     # Minimum number of observsations to run model
correct.bias = c(TRUE, FALSE)[1]                # Apply bias correction to occurance data 
r.dist = 100                                    # distance threshold for correcting road bias
bias.size = 0.10                                # subsampling size for correcting road bias
model.select = c(TRUE, FALSE)[2]                # Apply model selection
uncertainty = c(TRUE, FALSE)[2]                 # estimate spatial uncertainty
multi.thread = c(TRUE, FALSE)[2]                # parallel pcrocessing of predict
ncore = parallel::detectCores()-1               # number of cores for parallel processing
which.dat = c("all", "GPS", "opportunistic")[3] # what data to use for mammals
bird.rank = c(0:4)[1]                           # what weights to run for bird models

out.plots = "SDM_model.pdf"           # name of pdf plots output 
out.mdl = "SDM_model.RData"           # name of saved R model object ".RData"
out.val = "SDM_validation.txt"        # name of validation and parameter report
out.sdm = "spp_results.tif"           # name of output SDM probably and classification rasters 

#********************************************
# needed functions
source(file.path(root, "code", "occurrence.threshold.R")) 
source(file.path(root, "code", "accuracy.R"))
source(file.path(root, "code", "spatial.uncertainty.R")) 

pfun <- function(object, newdata, i = 2) {
  predict(object, data = newdata)$predictions[,i]
}

se.fun <- function(object, newdata, i = 2) {
  se <- predict(object, data = newdata, type = "se", se.method = "infjack")
  data.frame(p=se$predictions[,i], se=se$se[,i])  
}

dup.pairs <- function(x) {
  dx <- duplicated(x)
    other.dups <- which(dx)
      not.and.first <- which(!dx)
        ind.dups <- not.and.first[which(x[not.and.first] %in% x[other.dups])]
          xo <- x[other.dups]
  vc <- vector()
  for (cc in 1:length(ind.dups)) {
    vc <- c(vc,ind.dups[cc], other.dups[which(xo %in% x[ind.dups[cc]])])
  }
  return(vc)
}

rmse <- function(y,x) { sqrt(mean((y - x)^2)) }

brier.score <-  function(x, y) { mean((x - y)^2) }

#********************************************
# read raster covariates stack and create
# 30m reference raster for pseudo-absences

#**** 30m parameters
# parms <- c("climate_pca.tif", "cropland.tif", "developed.tif", "dist_dev.tif", 
#            "dist_forest.tif", "dist_lake.tif", "dist_roads.tif", "dist_stream.tif", 
#            "forest.tif", "forest_cover.tif", "forest_height.tif", "forest_type_pct.tif", 
#            "grassland.tif", "hli.tif", "LAI_trend.tif", "scosa.tif", "srr.tif", "tpi.tif", 
#            "tri.tif", "wetland.tif")        
parms <- c("climate_pca.tif", "crops.tif", "development_multiscale.tif", "dist_dev.tif", 
           "dist_forest.tif", "dist_lake.tif", "dist_roads.tif", "dist_stream.tif", 
		   "forest_cover.tif", "forest_heights.tif", "forest_pct.tif", 
		   "grasslands_multiscale.tif", "hli.tif", "LAI_trend.tif", "scosa.tif", "srr.tif", 
            "tpi.tif", "tri.tif", "water.tif", "xcoord.tif", "ycoord.tif")
r <- rast(file.path(dat.dir, parms))	   

bdy <- st_transform(st_read(file.path(dat.dir, "Croatia.gpkg"), "boundary"), prj)
ref <- rast(file.path(dat.dir, "mask100m.tif"))
dev <- rast(file.path(dat.dir, "developed_mask.tif"))
  dev <- mask(crop(dev, ref), ref)

pa.ref <- rast(ext(ref), resolution = 300, crs=crs(ref))
  pa.ref[] <- 1
    pa.ref <- mask(pa.ref, vect(bdy))
biogeo <- st_transform(st_cast(st_read(file.path(dat.dir, "Croatia.gpkg"), 
                       "biogeographic_regions"), "POLYGON"), prj)

# collapse "pannonian" into "continental"   
biogeo$short_name[which(biogeo$short_name == "pannonian")] <- "continental" 

if(is.null(nsample)){
  rn <- length(ref[!is.na(ref)][,1])
  nsample <- round(rn * sfract, 0)
}

#********************************************
#********************************************
# Model Loop
if(group == "birds"){
  spp.scores <- read.csv(file.path(mdl.dir, "bird_list.csv"))
    spp.scores$priority[which(is.na(spp.scores$priority))] <- 0
      mdl.spp <- spp.scores[which(spp.scores$priority == bird.rank),]$species
  spp.names <- data.frame(spp = st_layers(file.path(mdl.dir, spp.dat))$name,
                          cts = st_layers(file.path(mdl.dir, spp.dat))$features)
    spp.names <- spp.names[which(spp.names$cts >= min.ct),]$spp
      spp.names <- spp.names[which(spp.names %in% mdl.spp)]
} else {
  spp.names <- data.frame(spp = st_layers(file.path(mdl.dir, spp.dat))$name,
                          cts = st_layers(file.path(mdl.dir, spp.dat))$features)
    spp.names <- spp.names[which(spp.names$cts >= min.ct),]$spp
}
fin.mds <- basename(list.dirs(mdl.dir))[-1]
  fin.mds <- which(spp.names %in% fin.mds)
    if(length(fin.mds) > 0)
      spp.names <- spp.names[-fin.mds]

#*** Subset models for running multiple R sessions
# spp.names <- spp.names[1:65]
# spp.names <- spp.names[66:130]

ct=0    
  for(i in spp.names) {
  ct=ct+1
    cat("\n", "Processing", basename(i), "\n")
	
    #***************************************************
    # read and build data	
	spp <- st_transform(st_read(file.path(mdl.dir, spp.dat), i, quiet = TRUE), prj)
	  st_geometry(spp) <- "geometry"
        if(group == "mammals"){
		  cat("**** Using:", which.dat, "data for", i, "\n")
          if(which.dat != "all") {
            spp <- spp[spp$sample_type == which.dat,]
          }
        }
    didx <- which(duplicated(suppressWarnings(extract(ref, spp, cells=TRUE)$cell)))
      if(length(didx) > 0) {
	    message("removing ", length(didx), " duplicate cell locations")
	    spp <- spp[-didx,]
      }		
    if(nrow(spp) < min.ct) {
	  message("Not engough observations to run model, skipping to next")
  	  next
	}
    dn <- file.path(mdl.dir, i)
      dir.create(dn, showWarnings = FALSE)
        setwd(dn)

    report <- list()
	  report[["plots.out"]] <- paste0("plots written to ", file.path(i, out.plots))
	  report[["mdl.out"]] <- paste0("R Model written to ", file.path(i, out.mdl))
	  report[["rep.out"]] <- paste0("Validation report written to ", file.path(i, out.val))
      report[["train.out"]] <- paste0("training data written to ", file.path(i, out.sdm))
      report[["sdm.out"]] <- paste0("SDM rasters written to ", file.path(i, out.plots))
      report[["nboot"]] <- paste0("Number of Random Forests Bootstrap replicates: ", nboot)
	  report[["pasamp"]] <- paste0("Number of pseudo-absence observations: ", nsample)

    #***************************************************
    # occurance debias
    if(correct.bias){
	rd.wts <- c(0.001, 0.1, 0.05, 0.8)
	  cat("***************************************************", "\n")
      cat("**** Correcting for sampling bias", "\n")
      cat("***************************************************", "\n")
      rds <- rast(file.path(dat.dir, "dist_roads.tif")) 
        drds <- extract(rds, spp)[,-1]
          na.idx <- unique(which(is.na(drds), arr.ind = TRUE)[,1])
      	    if(length(na.idx) > 0 ) {
      	      spp <- spp[-na.idx,]
      	      drds <- drds[-na.idx,]
      	    }
      norg <- nrow(spp)			
      d.idx <- unique(which(drds$dist_hwy <= r.dist |
	                  drds$dist_state_rds <= r.dist |
                      drds$dist_county_rds <= r.dist |
      	              drds$dist_local_rds <= r.dist))	  
        if(length(d.idx) > 10) {	  
          spp.rds <- spp[d.idx,]
          spp <- spp[-d.idx,]
          drds <- drds / max(drds)
          wts <- apply(drds[d.idx,], MARGIN=1, FUN=\(x) {
                         x <- sum(x * rd.wts) })
            n = nrow(spp.rds)
	    	s = round(n * bias.size, 0)
	    	  if(s == 0) s = 1
            samp.idx <- sample(1:n, s, prob = wts/sum(wts))
			  spp.bias <- spp.rds[-samp.idx,] 
	    	  spp.rds <- spp.rds[samp.idx,]			
			# rds <- sum(rds, na.rm=TRUE)
			#   rds <- rds / global(rds, max, na.rm=TRUE)[,1]
			#plot(rds, maxcell=1000000, legend=FALSE, col=rev(grDevices::terrain.colors(20)), breaks = 20)
			#  plot(st_geometry(spp), pch=20,cex=0.8, add=TRUE)
            #    plot(st_geometry(spp.bias), pch=20,cex=1,col="red",add=TRUE)
			#	plot(st_geometry(spp.rds), pch=20,cex=1.5,col="green",add=TRUE)
	    	spp <- rbind(spp, spp.rds)
			report[["bias.corr"]] <- paste0("Removed ", norg-nrow(spp),  " bias points - ", round(nrow(spp)/norg*100,0), " percent")
			cat("**** Removed", norg-nrow(spp),  "bias points -", 100 - round(nrow(spp)/norg*100,0), "percent", "\n")
          }
	  }

	  cat("***************************************************", "\n")
      cat("**** Probabilistic Random Forests for:", i, "\n")
      cat("****   Model has ", nrow(spp), " prevleance observations", "\n")  
      cat("****   Running ", ct, " of ", length(spp.names), " models", "\n") 
      cat("***************************************************", "\n") 	  
		
    #***************************************************
    # mask to observed biogeographic regions
    spp.biogeo <- table(factor(suppressWarnings(sf::st_intersection(spp, biogeo)$short_name),
	                    levels=c("mediterranean", "alpine", "continental")))
    # spp.biogeo <- unique(suppressWarnings(sf::st_intersection(spp, biogeo)$short_name))
	  if(any(spp.biogeo < min.ct)) {
	    cat("***************************************************", "\n")
        cat("**** Subseting data to ", length(unique(biogeo$short_name)),  " biogeographic regions", "\n")
        cat("***************************************************", "\n")
        spp.biogeo <- names(which(spp.biogeo >= min.ct)) 		
		if(length(spp.biogeo) > 0) {
		biogeo.sub <- biogeo[which(biogeo$short_name %in% spp.biogeo),]
          spp.mask <- mask(ref, vect(biogeo.sub))
            spp.parms <- mask(r, spp.mask)
		      spp.ref <- mask(pa.ref, vect(biogeo.sub))
		} else {
          message("Not enough observations per geographic region, skipping model")
		    next
        }		
      } else {
	    cat("***************************************************", "\n")
        cat("**** Using full geographic extent of data", "\n")
        cat("***************************************************", "\n") 		  
        spp.mask <- ref 
		  spp.parms <- r
		spp.ref <- pa.ref
	  }
	  
    #***************************************************
    # Create pseudo-absence
	nsample <- round(nrow(spp) * 2.5, 0)
      #if(is.null(nsample)){
      #  rn <- length(ref[!is.na(ref)][,1])
      #  nsample <- round(rn * sfract, 0)
      #}
    pa <- pseudo.absence(spp, n = nsample, KDE=FALSE, ref = spp.ref, sigma=bw.method)
	  report[["kde.bandwidth"]] = pa$sigma
      report[["kde.plugin"]] = bw.method 
      pa$sample <- cbind(pa$sample, 
          data.frame(y=rep(0,nrow(pa$sample)),
          obs_type="pseudo abs", 
          species = unique(spp$species) ) )
        pa$sample <- pa$sample[,-1]
      spp$y <- 1
        spp$obs_type <- "occurance"
          spp <- spp[,c("y", "obs_type", "species")]
	st_crs(pa$sample) <- st_crs(spp)  
	spp <- rbind(spp, pa$sample)
      remove(pa)
	    gc()

    #***************************************************
    # Extract raster covariates
    x <- extract(spp.parms, spp, cells = TRUE)[,-1]
	  cells <- x$cell
	    x <- x[,-ncol(x)]
    didx <- which(spp$y == 0 & duplicated(cells))
      if(length(didx) > 0) {
	    message("removing ", length(didx), " duplicate cell locations")
	    spp <- spp[-didx,]
		x <- x[-didx,]
      }			
      na.idx <- unique(which(is.na(x), arr.ind = TRUE)[,1])
        if(length(na.idx) > 0) {
		  cat("removing", length(na.idx), "NA observsations", "\n")
	      x <- x[-na.idx,]
	      spp <- spp[-na.idx,]	
	    }

    #***************************************************
    # pairwise and multi collinearity screening 
    all.vars <- names(x)
      cl.vars <- try(collinear(x, p = 0.80, nonlinear = FALSE, 
	                 p.value = 0.001))
  	  if(exists("cl.vars")){
        if(length(cl.vars) > 0)
          x <- x[,-which(names(x) %in% cl.vars)]
  	  }	
      mc <- multi.collinear(x, p = 0.05, perm = TRUE, n = 99)
        mc.vars <- mc[which(mc$frequency > 10 ),]$variables 
          if(length(mc.vars) > 0){
            x <- x[,-which(names(x) %in% mc.vars)]
  		}	
      cat("The following variables were dropped:",
        all.vars[which(is.na(match(all.vars, names(x))))], "\n")
	  
	# Combine presence and pseudo-absence, write shapefile 
	spp <- cbind(spp, x)
	  dat <- st_drop_geometry(spp)[,-c(2:3)]
	    dat$y <- factor(dat$y)
	
	report[["dropped.vars"]] <- all.vars[which(is.na(match(all.vars, names(x))))]
      report[["pres.n"]] <- nrow(dat[dat$y == "1",])
        report[["null.n"]] <- nrow(dat[dat$y == "0",])

    #***************************************************
    # Model/parameter selection
	# Test variable significance p-value and remove negative importance
	rf.exp <- ranger(y ~ ., data = dat, probability = TRUE, 
                     num.trees = nboot, importance="impurity_corrected")
      ( imp.p <- as.data.frame(importance_pvalues(rf.exp, method = "altmann", 
                             formula = y ~ ., data = dat)) )
        nvars <- rownames(imp.p)[which(imp.p$pvalue < 0.05 & imp.p$importance > 0)] 
          dat <- data.frame(y=factor(dat$y), x[,nvars])
   non.sig <- rownames(imp.p)[which(imp.p$pvalue >= 0.05 & imp.p$importance >= 0)]  
      if(length(non.sig) == 0) non.sig = NULL
    if(model.select) {	  
      rf.sel <- rf.modelSel(xdata = dat[2:ncol(dat)], ydata = dat$y, imp.scale = "mir", 
                            r = c(0.10, 0.25, 0.50, 0.75, 0.90), seed = 42, ntree = nboot) 
        sel.vars <- rf.sel$selvars
          dat <- data.frame(y=dat$y, dat[,sel.vars])
    } else {
	  sel.vars <- names(dat)[-1]
	}
	  report[["nonsig.vars"]] <- non.sig 
    report[["sel.vars"]] <- sel.vars

    #***************************************************
    # Fit model		
    rf.fit <- ranger(y = dat$y, x = dat[,2:ncol(dat)], probability = TRUE, 
                     num.trees = nboot, importance="permutation",
				     replace = TRUE, keep.inbag=TRUE)     
      
    # Optimize fit on NULL class prediction variance
    p <- pfun(rf.fit, dat)
    pidx <- which(dat$y == "0" & p > 0.50)  
      if(length(pidx) > 0) {
	    dat <- dat[-pidx,]
        rf.fit <- ranger(y = dat$y, x = dat[,2:ncol(dat)], probability = TRUE, 
                         num.trees = nboot, importance="permutation",
						 replace = TRUE, keep.inbag = TRUE)     
      }

    ## Check Bootstrap sample sizes
    # inbag <- do.call(cbind, rf.fit$inbag.counts)
    #   unique(colSums(inbag[which(dat$y == "0"),]))
    #   unique(colSums(inbag[which(dat$y == "1"), ]))

    #***************************************************
    # Validation of fit and performance

	cat("***************************************************", "\n")
    cat("**** Cross-validation of fit and performance for:", i, "\n")
    cat("***************************************************", "\n") 

	se <- data.frame(y=as.numeric(as.character(dat$y)), se.fun(rf.fit, dat))
	
	#**** Confusion matrix fit
    #**** Global and local log loss, log likelihood fit and Brier score
	report[["validation"]] <- accuracy(table(se$y, ifelse(se$p > 0.6, 1, 0))) 
    report[["log.loss"]] <- logLoss(y = se$y, p = se$p) 		
    report[["Brier.score"]] <- brier.score(se$y, se$p)	
		
    #**** cross-validation w/ 10% Bootstrap
	# ncv = 99
	# cv <- data.frame(matrix(ncol = 6, nrow = ncv))
	#   names(cv) <- c("log.loss", "error", "brier", "sensitivity", "specificity", "kappa")
	# for(j in 1:ncv){
	#   samp.size <- round(length(which(dat$y %in% "1")) * 0.10,0)
	#     if(samp.size < 20) samp.size = 20
	#   idx1 <- sample(which(dat$y %in% "1"), samp.size) 
	#   idx0 <- sample(which(dat$y %in% "0"), samp.size) 
	#   withold <- rbind(dat[idx0,], dat[idx1,])
	#   p.idx <- which(withold$y == 1)
	#   dat.sub <- dat[-c(idx1, idx0),]
    #   rf.cv <- ranger(y = dat.sub$y, x = dat.sub[,2:ncol(dat)], probability = TRUE, 
    #                      num.trees = nboot, importance="permutation",
	# 					 replace = TRUE, keep.inbag = TRUE)   
	#     cv.se <- data.frame(y=as.numeric(as.character(withold$y)), 
	# 	                    suppressWarnings(se.fun(rf.cv, withold)))
	# 	cv.val <- accuracy(table(cv.se$y, ifelse(cv.se$p > 0.6, 1, 0)))
	# 
	# 	cv[j,] <- c(as.numeric(logLoss(y = cv.se$y, p = cv.se$p)),
	# 	            brier.score(cv.se$y, cv.se$p),
	# 	            rf.cv$prediction.error, 
	# 	            cv.val$sensitivity,
	# 	            cv.val$specificity,
	# 	            cv.val$kappa) 
	# }
	#           report[["cv.logloss"]] <- median(cv$log.loss)
	#         report[["cv.error"]] <- median(cv$error)
	#       report[["cv.brier"]] <- median(cv$brier)
	#     report[["cv.sensitivity"]] <- median(cv$sensitivity)
	#   report[["cv.specificity"]] <- median(cv$specificity)
	# report[["cv.kappa"]] <- median(cv$kappa)
	
    #***************************************************
    # Model prediction
	
	cat("***************************************************", "\n")
    cat("**** Making spatial predictions for:", i, "\n")
	cat("**** ", ct, " of ", length(spp.names), " models", "\n") 
    cat("***************************************************", "\n") 	  	
	
	# Threshold sensitivity test (mean of multiple metrics)
	pt <- seq(0.40, 0.80, 0.05)
    sum.ss <- suppressWarnings(occurrence.threshold(rf.fit,  p=pt,  
                               class = "1", type = "sum.ss")) 
    kappa.ss <- suppressWarnings(occurrence.threshold(rf.fit,  p=pt,  
                                 class = "1", type = "kappa"))
    logloss.ss <- suppressWarnings(occurrence.threshold(rf.fit,  p=pt,  
                                 class = "1", type = "logloss"))
     thresholds <- c(round(mean(c(sum.ss$prob.threshold, 
	                kappa.ss$prob.threshold,
                    logloss.ss$prob.threshold)),3),
					sum.ss$prob.threshold, kappa.ss$prob.threshold, 
					logloss.ss$prob.threshold) 
      names(thresholds) <- c("threshold", "sum.ss", "kappa", "log.loss")
    report[["p.threshold"]] <- thresholds 

    # sdm <- spp.parms[[1]]
	#   sdm[] <- NA
	# spp.parms <- spp.parms[[which(names(spp.parms) %in% sel.vars)]]   
    #   spp.parms <- as.data.frame(spp.parms, cells = TRUE, na.rm = TRUE)
	# sdm.pred <- predict(rf.fit, data=spp.parms[,-1], type="response")$predictions[,2]
	#   sdm[spp.parms$cell] <- sdm.pred     
    #
	#  rfp <- randomForest(y = dat$y, x = dat[,-1], ntree = nboot)
	#    if(multi.thread) {
	#      sdm <- predict(spp.parms, rfp,  type="prob", index=2, na.rm = TRUE, 
	#	                 cores = ncore, cpkgs = "randomForest")
    #    } else {
	#      sdm <- predict(spp.parms, rfp,  type="prob", index=2, na.rm = TRUE)	
	#    }
	spp.parms <- spp.parms[[which(names(spp.parms) %in% sel.vars)]] 
	  if(nlyr(spp.parms) < 3 ) {
	    message("Not enough selected paramters to support an estimate")
	    next
	  }
	 sdm <- predict(spp.parms, rf.fit, fun = pfun, na.rm = TRUE)
	   names(sdm) <- "SDM"
	     sdm <- mask(ifel(dev == 1, 0, sdm), sdm)
	 sdm.class <- ifel(sdm >= report$p.threshold[1], 1, 0)
	   names(sdm.class) <- "SDM_CLASS"
	     sdm <- c(sdm, sdm.class)	
      if(uncertainty) {
	    se.fun <- \(...) { predict(...)$se[,2] }
	    se <- predict(spp.parms, rf.fit, fun = se.fun, na.rm = TRUE, 
		              type = "se", se.method = "infjack")
          se <- c(se, y - (se * 1.96), y + (se * 1.96) )
            names(ci95) <- c("std.err", "lower.ci", "upper.ci")
		sdm <- c(sdm, se)
	  }  
    writeRaster(sdm, "sdm_estimates.tif", overwrite = TRUE)	  
      st_write(spp, "model_data.gpkg", "traning_data", append=FALSE)

    #***************************************************
    # Create model plots	
	nvars <- names(ranger::importance(rf.fit))[which(ranger::importance(rf.fit)>0)]
      imp <- ranger::importance(rf.fit)
        imp.names <- names(imp)
          idx <- order(imp)		
            imp <- imp[idx]
              imp.names <- imp.names[idx]
                imp = imp / max(imp)
			  
	pdf(out.plots, height=10, width=10)		
      # Plot importance  
	  dotchart(imp, imp.names, pch=19, main="Variable Importance")

      # Plot probabilities
	  pal <- ifelse(spp$y == "1", "red", "blue")      
      plot(sdm[[1]], main=paste0("Estimated Probabilities for ", unique(spp$species)))	  
	    plot(spp["y"], pch=20, cex=0.75, col=pal, add=TRUE)
	      legend(4620000, 2250000,legend=c("presence", "absence"),
		         pch=c(20,20), col=c("red", "blue"))   
      plot(sdm[[2]], main=paste0("Classified probs for ", unique(spp$species)))	  
	    plot(spp["y"], pch=20, cex=0.75, col=pal, add=TRUE)
	      legend(4620000, 2250000,legend=c("presence", "absence"),
		         pch=c(20,20), col=c("red", "blue"))   
	  
	#**** validation plots
	# Cross-validation log-loss and error  
    #  suppressWarnings(plot(sort(cv.logloss),type="line",
	#                   main="Cross-validation logloss",
	#				   ylab="Log Loss", xlab=""))
    #                   abline(h=gll, lty=3)
    #  suppressWarnings(plot(sort(cv.error),type="line",
	#	               main="Cross-validation error",
	#				   ylab="Error", xlab=""))
    #                   abline(h=rf.fit$prediction.error , lty=3)
	  
	  # Plot standard errors
      pal <- ifelse(dat$y == "1", "blue", "red")
      classes <- ifelse(se$p > 0.5, "present", "absent")
      plot(se$p, se$se, col = pal, main="Standard Errors", 
           pch = 19, xlab = "Predicted probability",
           ylab = "Standard error")
      
	  # Plot confidence interval             
      scalar <- (1.96 * se$se ) 
      ci95 <- data.frame(
      p = se[,1], std.err=se[,2],
      	lower.ci = se[,1] - scalar, 
          upper.ci = se[,1] + scalar)
      matclr <- ifelse(dat$y == 1, "blue", "red")
      sort.idx <- sort(ci95$p, index.return = TRUE)$ix	  
      plot(ci95$p[sort.idx], type="p", xaxt="n", col=matclr[sort.idx],
           pch=20, cex=0.5, ylab="probabilities", xlab="observations (n=1277)", 
      	 main="probability estimate with 95% confidence intervals",  
           ylim=c(min(ci95[,3:4]), max(ci95[,3:4])))
        lines(ci95$upper.ci[sort.idx], col="grey", lty=2)
        lines(ci95$lower.ci[sort.idx], col="grey", lty=2)
      legend("bottomright", legend=c("Presence", "Absence"),
             pch=20, col=c("blue", "red")) 

	  # Plot Log Loss 
      ll <- logLoss(y = se$y, p = se$p, global=FALSE) 
        lld <- density(ll$log.loss)
          plot(lld, type="n", main="prediction log loss", xlab="Log Loss", ylab="pdf")
            polygon(lld, col="blue")

      # Partial dependency plots (functional relationships)
      for(n in names(rf.fit$variable.importance)) {
        cat("Calculating partial plot for parameter:", n, "\n")  
        pp <- pdp::partial(rf.fit, pred.var = n, plot = TRUE, which.class = 1, 
      	                 prob = TRUE, train = dat[,-1], 
      			         levelplot = TRUE,chull = TRUE, quantiles = TRUE, 
      					 smooth=FALSE, plot.engine = "ggplot2")
        print(pp + ggtitle(paste0("Partial dependency for ", n)))
      }
	
      # Plot Threshold sensitivity test	
      par(mfrow=c(2,2))
        plot(sum.ss)
        plot(kappa.ss)  
	    plot(logloss.ss)  
    dev.off()
	
	# SDM report with used parameters and validation
	
	file.create(file.path(getwd(), out.val))
    write(paste0("Output files"), file=(file.path(getwd(), out.val)), append=TRUE)		
	  write(report$rep.out, file=(file.path(getwd(), out.val)), append=TRUE)
	  write(report$sdm.out, file=(file.path(getwd(), out.val)), append=TRUE)
	  write(report$train.out, file=(file.path(getwd(), out.val)), append=TRUE)
      write(report$plots.out, file=(file.path(getwd(), out.val)), append=TRUE)
      write(report$mdl.out, file=(file.path(getwd(), out.val)), append=TRUE)
	  write(paste0(), file=(file.path(getwd(), out.val)), append=TRUE)
      write(paste0("Model parameters"), file=(file.path(getwd(), out.val)), append=TRUE)	  
	  write(report$nboot, file=(file.path(getwd(), out.val)), append=TRUE)	  	  
	  write(report$nsamp, file=(file.path(getwd(), out.val)), append=TRUE)	  	  
	  write(paste0("Pseudo-absence KDE bandwidth: ", report$kde.bandwidth),
            file=(file.path(getwd(), out.val)), append=TRUE)	  
	  write(paste0("Pseudo-absence KDE method: ", report$kde.plugin),
            file=(file.path(getwd(), out.val)), append=TRUE)	  
	  write(paste0("Number of presence obs: ", report$pres.n),
            file=(file.path(getwd(), out.val)), append=TRUE)	  
	  write(paste0("Number of Pseudo-absence obs: ", report$null.n),
            file=(file.path(getwd(), out.val)), append=TRUE)	  
	  write(paste0(), file=(file.path(getwd(), out.val)), append=TRUE)
      write(paste0("Model selection"), file=(file.path(getwd(), out.val)), append=TRUE)	  
      write(paste0("Dropped correlated variables: ", report$dropped.vars),
	        file=(file.path(getwd(), out.val)), append=TRUE)	
      write(paste0("Dropped non-significant variables: ", report$nonsig.vars),
	        file=(file.path(getwd(), out.val)), append=TRUE)	
      write(paste0("Selected parameters: ", report$sel.vars),
		    file=(file.path(getwd(), out.val)), append=TRUE)		 
	  write(paste0(), file=(file.path(getwd(), out.val)), append=TRUE)
      write(paste0("Model validation - fit"), file=(file.path(getwd(), out.val)), append=TRUE)
      write(paste0("Percent Correctly Classified: ", report$validation$PCC),
            file=(file.path(getwd(), out.val)), append=TRUE) 	  
      write(paste0("Kappa: ", report$validation$kappa),
            file=(file.path(getwd(), out.val)), append=TRUE) 	  
      write(paste0("Log loss: ", report$log.loss),
            file=(file.path(getwd(), out.val)), append=TRUE) 	  
      write(paste0("Brier score: ", report$Brier.score),
            file=(file.path(getwd(), out.val)), append=TRUE) 	  			
      write(paste0("Sensitivity: ", report$validation$sensitivity),
            file=(file.path(getwd(), out.val)), append=TRUE) 	  
      write(paste0("Specificity: ", report$validation$specificity),
            file=(file.path(getwd(), out.val)), append=TRUE) 	  
      write(paste0("Area Under ROC (AUC): ", report$validation$auc),
            file=(file.path(getwd(), out.val)), append=TRUE) 	  
      write(paste0("Type I error: ", report$validation$typeI.error),
            file=(file.path(getwd(), out.val)), append=TRUE) 	
      write(paste0("Type II error: ", report$validation$typeII.error),
            file=(file.path(getwd(), out.val)), append=TRUE) 	  		
      write(paste0("Matthew's correlation): ", report$validation$matthewsc),
            file=(file.path(getwd(), out.val)), append=TRUE) 	  
      write(paste0("Information Gain (signal to noise)): ", report$validation$gain),
            file=(file.path(getwd(), out.val)), append=TRUE) 	  
	  write(paste0(), file=(file.path(getwd(), out.val)), append=TRUE)
      # write(paste0("Model validation - performance (cross-validation)"), 
	  #       file=(file.path(getwd(), out.val)), append=TRUE)
      # write(paste0("Cross-validation error: ", report$cv.error),
	  #       file=(file.path(getwd(), out.val)), append=TRUE)	
      # write(paste0("Cross-validation log loss: ", report$cv.logloss),
	  #       file=(file.path(getwd(), out.val)), append=TRUE)	
      # write(paste0("Cross-validation Brier Score: ", report$cv.brier),
	  #        file=(file.path(getwd(), out.val)), append=TRUE)	
      # write(paste0("Cross-validation sensitivity: ", report$cv.sensitivity),
	  #        file=(file.path(getwd(), out.val)), append=TRUE)	
      # write(paste0("Cross-validation specificity: ", report$cv.specificity),
	  #        file=(file.path(getwd(), out.val)), append=TRUE)	
      # write(paste0("Cross-validation kappa: ", report$cv.kappa),
	  #        file=(file.path(getwd(), out.val)), append=TRUE)				
	  write(paste0(), file=(file.path(getwd(), out.val)), append=TRUE)
      write(paste0("Probability threshold for binomial classification"), 
	        file=(file.path(getwd(), out.val)), append=TRUE)
      write(paste0("Probability classification threshold ", 
	        c("mean: ", "sum.ss: ", "kappa: ", "log loss: "), 
			report$p.threshold), 
	        file=(file.path(getwd(), out.val)), append=TRUE)	  
	# file.show(file.path(getwd(), out.val))
	
    save.image(file.path(getwd(), out.mdl))
      # terra::tmpFiles(current=TRUE, orphan=TRUE, old=TRUE, remove=TRUE)
    remove(dat, rf.fit, prf, spp, spp.parms, sdm, rf.exp, rf.sel)   
  gc()	
  } # End model Loop

#********************************************
#********************************************
