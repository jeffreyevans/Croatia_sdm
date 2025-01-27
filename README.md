# Croatia SDM's version 1.0 (01/27/2025)

Croatia's Species Distribution Models. 

Description of model code files 

1. create_sdm_covariates.R - 

2. Download_GLASS_LAI.R - Downloads and processes Leaf Area Index from GLASS, returns trend covariate

3. process_telemetry_data.R - Creates a density-based sample of the telemetry data using an Autocorrelated Kernel Density Estimate 

4. croatia_sdm.R and croatia_telemetry_sdm.R - Species Distribution Model. The telemetry version creates a balanced subsample between telemetry and occurance data   

Other code

1. accuracy.R – Modification of rfUtilities function which provides validation metrics based on confusion matrix

2. check.packages.R - Checks if required libraries are installed and adds them if not in library. Also defines library environment and mirror. Can be used to simply add required libraries to current R environmnet. 

3. occurrence.threshold.R - Modification of rfUtilities function to include log loss as option

4. RangeOfVariability.R - Evaluates the range of variability (distributions) of covariates within bioregions

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
