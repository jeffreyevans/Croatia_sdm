# set site library
.Library.site <- file.path(chartr("\\", "/", R.home()), "library")
.libPaths(file.path(chartr("\\", "/", R.home()), "library"))

# set a CRAN mirror
local({r <- getOption("repos")
       r["CRAN"] <- "https://ftp.osuosl.org/pub/cran/"
       options(repos=r)})

# Check libraries, if not installed than add, else add to environment
chk.pkg <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg)) 
        install.packages(new.pkg, dependencies = TRUE)
    sapply(pkg, require, character.only = TRUE)
}

pkg <- c("sf", "spatialEco", "terra", "ranger", "ggplot2",  
         "randomForest", "rfUtilities", "pdp", "readxl", 
		 "geodata")
chk.pkg(pkg)
