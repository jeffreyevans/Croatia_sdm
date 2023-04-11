invisible(lapply(c("sf", "spatialEco", "terra", "ranger", 
          "ggplot2", "randomForest", "rfUtilities", "pdp"), 
		  require, character.only=TRUE))

set.seed(42)                          # Random seed, for reproducibility 

root = "C:/evans/Croatia"             # Main root directory 
setwd(root)                           # Set working directory to main root directory
  dat.dir <- file.path(root, "data")  # directory with raster covariates
  mdl.dir <- file.path(root, "SDM")   # root directory of model directories

# number of pseudo-absence samples, if nsample = NULL  
# then the defined sample fraction (sfract) will be used, 
# based on number 0f non-NA raster, to draw sample size.
# At 100m are 237,067 no-NA cell values
nsample = 1000  # or NULL
sfract = 0.01   # eg., 1% sample of raster cells (n=2371)                       

nboot = 1001                          # number of bootstrap replicates
bw.method = c("Diggle", "Scott")[1]   # Bandwidth plug-in (1st or 2nd order)

out.plots = "SDM_model.pdf"           # name of pdf plots output 
out.mdl = "SDM_model.RData"           # name of saved R model object ".RData"
out.val = "SDM_validation.txt"        # name of validation and parameter report
out.sdm = "spp_results.tif"           # name of output SDM probably and classification rasters 

#********************************************
# needed functions
source(file.path(root, "code", "occurrence.threshold.R")) 
source(file.path(root, "code", "accuracy.R")) 

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

#********************************************
# read raster covariates stack and create
# 500m reference raster for pseudo-absences

r <- rast(list.files(dat.dir, pattern = "tif$", 
          full.names = TRUE))

bdy <- st_read(file.path(dat.dir, "bdy.shp"))
ref <- rast(ext(bdy), resolution=500, crs=crs(r))
  ref[] <- rep(1,ncell(ref))
    ref <- mask(ref, vect(bdy))

if(is.null(nsample)){
  rn <- length(ref[!is.na(ref)][,1])
  nsample <- round(rn * sfract, 0)
}

#********************************************
#********************************************
# Model Loop

d <- list.dirs(mdl.dir)[-c(1)]
  for(i in d) {
    cat("\n", "Processing", basename(i), "\n")
    setwd(i)
      report <- list()
	  report[["plots.out"]] <- paste0("plots written to ", file.path(i, out.plots))
	  report[["mdl.out"]] <- paste0("R Model written to ", file.path(i, out.mdl))
	  report[["rep.out"]] <- paste0("Validation report written to ", file.path(i, out.val))
      report[["train.out"]] <- paste0("training data written to ", file.path(i, out.sdm))
      report[["sdm.out"]] <- paste0("SDM rasters written to ", file.path(i, out.plots))
      report[["nboot"]] <- paste0("Number of Random Forests Bootstrap replicates: ", nboot)
	  report[["pasamp"]] <- paste0("Number of pseudo-absence observations: ", nsample)
 
	spp <- st_read(list.files(i, "shp$")[1], quiet = TRUE)
      pres.idx <- grep("true abs", spp$obs_type) 
        if(length(pres.idx) > 0) spp <- spp[-pres.idx,]

    #********************
    # Create pseudo-absence
    pa <- pseudo.absence(spp, n=nsample, KDE=FALSE, ref = ref, 
                         sigma=bw.method)
	
	report[["kde.bandwidth"]] = pa$sigma
    report[["kde.plugin"]] = bw.method 
	  
      pa$sample <- cbind(pa$sample, 
          data.frame(y=rep(0,nrow(pa$sample)),
          obs_type="pseudo abs", 
          species=unique(spp$species),
          nobs=0,
          dttm=as.Date(Sys.time(), "%Y-%m-%d", tz = "CET")))
        pa$sample <- pa$sample[,-1]   
	st_crs(pa$sample) <- st_crs(spp) 
	spp <- rbind(spp, pa$sample)
      remove(pa, pres.idx)
	    gc()

    #********************
    # Extract raster covariates
    x <- extract(r, spp, cells = TRUE)[,-1]
	  cells <- x$cell
	    x <- x[,-ncol(x)]
      na.idx <- unique(which(is.na(x), arr.ind = TRUE)[,1])
        if(length(na.idx) > 0) {
		  cat("removing", length(na.idx), "NA observsations", "\n")
	      x <- x[-na.idx,]
	      spp <- spp[-na.idx,]	
	    }

    #********************
    # Remove raster cell duplicates, retaining 1 when present 		
    dup.cell.idx <- dup.pairs(cells)
      if(length(dup.cell.idx) > 0) {
        dup <- spp[dup.cell.idx,]
          dup$cell <- cells[dup.cell.idx]
		  dup$idx <- dup.cell.idx  
		rm.idx <- vector()
        for(u in unique(dup$cell)) {
		  dup.sub <- dup[dup$cell == u,]
		    if(length(unique(dup.sub$y)) > 1){		      
			  rm.idx <- append(rm.idx, dup.sub[dup.sub$y != "1",]$idx )
			} else {
			  rm.idx <- append(rm.idx, dup.sub[-1,]$idx )			  
			}
		}
		if(length(rm.idx) > 0){
		  spp <- spp[-rm.idx,]
		  x <- x[-rm.idx,]
		}  
		remove(dup, cells, dup.sub, rm.idx, u)
		  gc()
	  }	

    #********************
    # pairwise and multi collinearity screening 
    all.vars <- names(x)    
    cl.vars <- try(collinear(x, p = 0.85, nonlinear = FALSE, p.value = 0.001))
	  if(exists("cl.vars")){
        if(length(cl.vars) > 0)
          x <- x[,-which(names(x) %in% cl.vars)]
	  }	
    mc <- multi.collinear(x, p = 0.05, perm = TRUE, n = 99)
      mc.vars <- mc[which( mc$frequency > 5 ),]$variables 
        if(length(mc.vars) > 0) 
          x <- x[,-which(names(x) %in% mc.vars)]
    cat("The following variables were dropped:",
      all.vars[which(is.na(match(all.vars, names(x))))], "\n")
		
	# Combine presence and pseudo-absence, write shapefile 
	spp <- cbind(spp, x)
	  #st_write(spp, "pa_spp.shp", delete_layer = TRUE)
	dat <- st_drop_geometry(spp)[,-c(2:5)]
	  dat$y <- as.factor(dat$y)
	
	report[["dropped.vars"]] <- all.vars[which(is.na(match(all.vars, names(x))))]
      report[["pres.n"]] <- nrow(dat[dat$y == "1",])
        report[["null.n"]] <- nrow(dat[dat$y == "0",])

    #********************
    # Model/parameter selection
	# Test variable significance p-value and remove negative importance
	rf.exp <- ranger(y ~ ., data = dat, probability = FALSE, 
                     num.trees = nboot, importance="impurity_corrected",
					 sample.fraction = c(0.5, 0.5), replace = TRUE)
      imp.p <- as.data.frame(importance_pvalues(rf.exp, method = "altmann", 
                             formula = y ~ ., data = dat))
        nvars <- rownames(imp.p)[which(imp.p$pvalue < 0.05 & imp.p$importance > 0)] 
          dat <- data.frame(y=factor(dat$y), x[,nvars])

    rf.sel <- rf.modelSel(xdata = dat[2:ncol(dat)], ydata = dat$y, imp.scale = "mir", 
                       r = c(0.10, 0.25, 0.50, 0.75, 0.90), seed = 42, ntree = nboot) 
      sel.vars <- rf.sel$selvars
        dat <- data.frame(y=dat$y, dat[,sel.vars])

    non.sig <- rownames(imp.p)[which(imp.p$pvalue >= 0.05 & imp.p$importance >= 0)]  
      if(length(non.sig) == 0) non.sig = NULL

	report[["nonsig.vars"]] <- non.sig 
    report[["sel.vars"]] <- sel.vars

    #********************
    # Fit model		
    rf.fit <- ranger(y = dat$y, x = dat[,2:ncol(dat)], probability = TRUE, 
                     num.trees = nboot, importance="permutation",
				     sample.fraction = c(0.5, 0.5), replace = TRUE,
      				 keep.inbag=TRUE)     
      
    # Optimize fit on prediction variance
    p <- pfun(rf.fit, dat)
    pidx <- which(dat$y == "0" & p > 0.50)  
      if(length(pidx) > 0) {
	    dat <- dat[-pidx,]
        rf.fit <- ranger(y = dat$y, x = dat[,2:ncol(dat)], probability = TRUE, 
                         num.trees = nboot, importance="permutation",
						 sample.fraction = c(0.5, 0.5), replace = TRUE,
      				     keep.inbag = TRUE)     
      }

    ## Check Bootstrap sample sizes
    # inbag <- do.call(cbind, rf.fit$inbag.counts)
    #   unique(colSums(inbag[which(dat$y == "0"),]))
    #   unique(colSums(inbag[which(dat$y == "1"), ]))

    #********************
    #* Validation of fit and performance
	se <- data.frame(y=as.numeric(as.character(dat$y)), se.fun(rf.fit, dat))
	
	#**** Confusion matrix fit
	report[["validation"]] <- accuracy(table(se$y, ifelse(se$p > 0.6, 1, 0))) 
    report[["log.loss"]] <- logLoss(y = se$y, p = se$p) 		
		
    #**** Global and local log loss, log likelihood fit
    gll <- logLoss(y = se$y, p = se$p)    
    ll <- logLoss(y = se$y, p = se$p, global=FALSE) 

    #**** cross-validation w/ 10% Bootstrap
	cv.logloss <- vector()
	cv.error <- vector()
	for(j in 1:99){
	  idx1 <- sample(which(dat$y %in% "1"), round(length(which(dat$y %in% "1")) * 0.10,0)) 
	  idx0 <- sample(which(dat$y %in% "0"), round(length(which(dat$y %in% "1")) * 0.10,0)) 
	  withold <- rbind(dat[idx0,], dat[idx1,])
	  dat.sub <- dat[-c(idx1, idx0),]
      rf.cv <- ranger(y = dat.sub$y, x = dat.sub[,2:ncol(dat)], probability = TRUE, 
                         num.trees = nboot, importance="permutation",
						 sample.fraction = c(0.5, 0.5), replace = TRUE,
      				     keep.inbag = TRUE)     
	    cv.se <- data.frame(y=as.numeric(as.character(withold$y)), 
		                    se.fun(rf.cv, withold))
		cv.logloss[j] <- logLoss(y = cv.se$y, p = cv.se$p)
		cv.error[j] <- rf.cv$prediction.error 
	}
	report[["cv.logloss"]] <- median(cv.logloss)
	report[["cv.error"]] <- median(cv.error)
	
    #********************
    # Spatial prediction
	
	# Threshold sensitivity test (mean of multiple metrics)
	pt <- seq(0.40, 0.80, 0.05)
    sum.ss <- suppressWarnings(occurrence.threshold(rf.fit,  p=pt,  
                               class = "1", type = "sum.ss")) 
    kappa.ss <- suppressWarnings(occurrence.threshold(rf.fit,  p=pt,  
                                 class = "1", type = "kappa"))
    logloss.ss <- suppressWarnings(occurrence.threshold(rf.fit,  p=pt,  
                                 class = "1", type = "logloss"))
    report[["p.threshold"]] <- round(mean(c(sum.ss$prob.threshold, 
	                            kappa.ss$prob.threshold,
                                logloss.ss$prob.threshold)),3)
    
	sdm <- predict(r, rf.fit, fun=pfun, na.rm=TRUE)
	  names(sdm) <- "SDM"
    sdm.class <- ifel(sdm >= report[["p.threshold"]], 1, 0)
	  names(sdm.class) <- "SDM_CLASS"
	sdm <- c(sdm, sdm.class)  
    writeRaster(c(sdm, sdm.class), "spp_results.tif",
                overwrite = TRUE)	  
    st_write(spp, "traning_data.shp", delete_layer = TRUE)
	
    #********************
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
      suppressWarnings(plot(sort(cv.logloss),type="line",
	                   main="Cross-validation logloss",
					   ylab="Log Loss", xlab=""))
                       abline(h=gll, lty=3)
      suppressWarnings(plot(sort(cv.error),type="line",
		               main="Cross-validation error",
					   ylab="Error", xlab=""))
                       abline(h=rf.fit$prediction.error , lty=3)
	  
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
	
	file.create(file.path(i, out.val))
    write(paste0("Output files"), file=(file.path(i, out.val)), append=TRUE)		
	  write(report$rep.out, file=(file.path(i, out.val)), append=TRUE)
	  write(report$sdm.out, file=(file.path(i, out.val)), append=TRUE)
	  write(report$train.out, file=(file.path(i, out.val)), append=TRUE)
      write(report$plots.out, file=(file.path(i, out.val)), append=TRUE)
      write(report$mdl.out, file=(file.path(i, out.val)), append=TRUE)
	write(paste0(), file=(file.path(i, out.val)), append=TRUE)
    write(paste0("Model parameters"), file=(file.path(i, out.val)), append=TRUE)	  
	  write(report$nboot, file=(file.path(i, out.val)), append=TRUE)	  	  
	  write(report$nsamp, file=(file.path(i, out.val)), append=TRUE)	  	  
	  write(paste0("Pseudo-absence KDE bandwidth: ", report$kde.bandwidth),
            file=(file.path(i, out.val)), append=TRUE)	  
	  write(paste0("Pseudo-absence KDE method: ", report$kde.plugin),
            file=(file.path(i, out.val)), append=TRUE)	  
	  write(paste0("Number of presence obs: ", report$pres.n),
            file=(file.path(i, out.val)), append=TRUE)	  
	  write(paste0("Number of Pseudo-absence obs: ", report$null.n),
            file=(file.path(i, out.val)), append=TRUE)	  
	write(paste0(), file=(file.path(i, out.val)), append=TRUE)
    write(paste0("Model fit"), file=(file.path(i, out.val)), append=TRUE)	  
      write(paste0("Dropped correlated variables: ", report$dropped.vars),
	        file=(file.path(i, out.val)), append=TRUE)	
      write(paste0("Dropped non-significant variables: ", report$nonsig.vars),
	        file=(file.path(i, out.val)), append=TRUE)	
      write(paste0("Selected parameters: ", report$sel.vars),
		    file=(file.path(i, out.val)), append=TRUE)		 
	write(paste0(), file=(file.path(i, out.val)), append=TRUE)
    write(paste0("Model validation - fit"), file=(file.path(i, out.val)), append=TRUE)
      write(paste0("Percent Correctly Classified: ", report$validation$PCC),
            file=(file.path(i, out.val)), append=TRUE) 	  
      write(paste0("Kappa: ", report$validation$kappa),
            file=(file.path(i, out.val)), append=TRUE) 	  
      write(paste0("Log loss: ", report$log.loss),
            file=(file.path(i, out.val)), append=TRUE) 	  
      write(paste0("Sensitivity: ", report$validation$sensitivity),
            file=(file.path(i, out.val)), append=TRUE) 	  
      write(paste0("Specificity: ", report$validation$specificity),
            file=(file.path(i, out.val)), append=TRUE) 	  
      write(paste0("Area Under ROC (AUC): ", report$validation$auc),
            file=(file.path(i, out.val)), append=TRUE) 	  
      write(paste0("Type I error: ", report$validation$typeI.error),
            file=(file.path(i, out.val)), append=TRUE) 	
      write(paste0("Type II error: ", report$validation$typeII.error),
            file=(file.path(i, out.val)), append=TRUE) 	  		
      write(paste0("Matthew's correlation): ", report$validation$matthewsc),
            file=(file.path(i, out.val)), append=TRUE) 	  
      write(paste0("Information Gain (signal to noise)): ", report$validation$gain),
            file=(file.path(i, out.val)), append=TRUE) 	  
	write(paste0(), file=(file.path(i, out.val)), append=TRUE)
    write(paste0("Model validation - performance (cross-validation)"), 
	      file=(file.path(i, out.val)), append=TRUE)
      write(paste0("Cross-validation error: ", report$cv.error),
	        file=(file.path(i, out.val)), append=TRUE)	
      write(paste0("Cross-validation log loss: ", report$cv.logloss),
	        file=(file.path(i, out.val)), append=TRUE)	
	write(paste0(), file=(file.path(i, out.val)), append=TRUE)
      write(paste0("Probability classification threshold: ", report$p.threshold), 
	        file=(file.path(i, out.val)), append=TRUE)	  
	# file.show(file.path(i, out.val))
	
    save.image(file.path(i, out.mdl))
  } # End model Loop

#********************************************
#********************************************
