# Croatia SDM's version 1.0 (01/27/2025)

Croatia's Species Distribution Models. 

Description of model code files 

1. create_sdm_covariates.R - Creates spatial covariates for SDM 

2. process_telemetry_data.R - Creates a density-based sample of the telemetry data using an Autocorrelated Kernel Density Estimate 

3. croatia_sdm.R and croatia_telemetry_sdm.R - Species Distribution Model. The telemetry version creates a balanced subsample between telemetry and occurance data   

4. process_timeseries.R - Applies a temporal decomposition model for removing seasional periodicity

5. timseries_correlation.R - Derives temporal trend for LAI (Kendall Tau statistic)

6. occurrence.threshold.R - Modification of rfUtilities function to include log loss as option in deriving prevliance probability threshold for creating binomial presence/absence rasters

7. accuracy.R – Modification of rfUtilities function which provides validation metrics based on confusion matrix

8. spatial.uncertainty.R - Applies and Infinitesimal Jackknife (Wager et al., 2014) to calculate standard errors, lower and upper 95% confidence interval(s). Can output spatial uncertanity rasters as well. 

Utility code

1. check.packages.R - Checks if required libraries are installed and adds them if not in library. Also defines library environment and mirror. Can be used to simply add required libraries to current R environmnet. 

2. Download_GLASS_LAI.R - Downloads and processes Leaf Area Index from GLASS, returns trend covariate

3. hli.R - Modification of Heat Load Index model

4. RangeOfVariability.R - Evaluates the range of variability (distributions) of covariates within bioregions

5. zip_by_rank.R - Creates zip file(s) by species rank

Contact:

Mate Zec
Southeast Europe Renewable Energy Siting Specialist
The Nature Conservancy
c/o: TNC Brussels office
Av. des Arts 44
1040 Brussels
Belgium
mate.zec@tnc.org

Jeffrey S. Evans, Ph.D.
Senior Landscape Ecologist & Biometrician 
The Nature Conservancy | Global Protect, Science 
jeffrey_evans@tnc.org
