# Discretized US Drought Data to Support Statistical Modeling

Summary: Drought is a costly and disruptive natural disaster, with widespread implications for agriculture, wildfire, and urban planning.  We present a novel data set on US drought built to enable computationally efficient spatio-temporal statistical and probabilistic models of drought. We converted drought data obtained from the widely-used US Drought Monitor (USDM) from continuous shape files to a 0.5 degree regular lattice. These data cover the Continental US from 2000 to mid-2022. Known environmental drivers of drought include those obtained from the North American Land Data Assimilation System (NLDAS-2), US Geological Survey (USGS) streamflow data, and National Oceanic and Atmospheric Administration (NOAA) teleconnections data. The USGS streamflow data is itself a new gridded data product, aggregating point-referenced stream discharges from across the US to a common lattice using watersheds to combine nearby stream data. The resulting data set permits statistical and probabilistic modeling of drought with explicit spatial and/or temporal dependence.  Such models could be used to forecast short-range and even season-to-season future droughts with uncertainty, extending the reach and value of the current US Drought Outlook produced by the National Weather Service Climate Prediction Center. 

The complete data set can be found at https://doi.org/10.5061/dryad.g1jwstqw7 . Each row in this data set corresponds to a grid location in a given week. The columns are the different drought and environmental variables. The variables were processed and generated using the scripts contained in this repository which are described below. 

This repository contains R and Matlab scripts that were used to generate a discretized data set on United State drought status with other relevant environmental variables.

RasterizeUSDMforDryad.R is a R script that imports shapefiles of drought from https://droughtmonitor.unl.edu/DmData/GISData.aspx, and rasterizes them to a 0.5 degree grid and saves a flat file csv.

TeleconnectionsForDryad.R is an R script that imports raw data on four teleconnections (PNA, ENSO, NAO, AO) from cpc.ncep.noaa.gov and converts them into weekly data and exports a flat file csv.

create_NLDAS_files_weekly_05deg.m is a Matlab script that imports netCDF data files from https://disc.gsfc.nasa.gov/datasets?keywords=NLDAS and rasterizes the variables to a weekly 0.5 degree grid.

The variables related to streamflow were generated using the following scripts:

1.  Data Download & Preprocessing

    a.	Downloaded text files with daily streamflow data (in cubic meters per second) from USGS website (https://waterdata.usgs.gov/nwis/dv/?referred_module=sw) for each hydrologic region. 
    
    b.	From the same website, downloaded the site geographic information for each hydrologic region as text files. 
    
    c.	Used Unix text editor to remove header information from daily flow data for each region, using ‘batchClean.sh’ on directory with daily flow text files.
    
    d.	Ran ‘readLatLonFindGrid.m’ to match gauge site geographic information to the model grids that are used in our statistical model, or the ‘common lattice’. This script creates a file called gridGaugeAlign.mat.
    
2.	readStreamflowDataFromTxt.m 
This script reads in the daily flow text files downloaded from USGS (with headers removed) and calculates the average 7-day, 14-day and 28-day flows and flow percentiles for each region. These flows are stored in MATLAB data structures. Each data structure contains an array for each site ID that is available within the region. This script requires gridGaugeAlign.mat data to run. The script can also calculate the percentage of missing data.

3.	getStreamflowModelGrid.m 
This script iterates through all regional streamflow data MATLAB files (created from readStreamflowDataFromTxt.m) and matches the site IDs with the model grids and reorganizes the data. Two MATLAB data structures are created for each temporal aggregation (7-day, 14-day, and 28-day). One contains a table for each model grid and stores the streamflow data for each individual gauge and the other contains a vector array for each model grid that contains the mean flow from all available gauges within the grid. This script also records the number of datapoints used in each calculation, from 2000 to 2021. This script requires gridGaugeAlign.mat and modelGridInfo.mat data files to run.

4.	exportData2csv.m
This script reads in the mean flow MATLAB data structures produced by getStreamflowModelGrid.m and creates a table that can be ingested by the statistical model. Each year of tabular data is exported to a csv file, called flowData(YEAR).csv. A value of NaN is used when there is no gauge data.
5.	fillMissingFlows.m 
This script sorts through missing flow values from the csv files produced by exportData2csv.m for each model grid and each date and calculates the number of gauges within the HUC-8/6/4 watersheds that could be used to gap-fill. It also calculates the distance between the centroid of the model grid and the gauge locations, in decimal degrees. This script depends on the data structures produced by getStreamflowModelGrid.m. The fill value information is stored in a MATLAB data structure for each year and each level of temporal aggregation. The script requires ‘gridGaugeAlign.mat’, ‘modelGridInfo.mat’, ‘stationHUCs.mat’ and ‘gaugeLocations.csv’ to run.
6.	getFillValsCompile.m 
This script uses the data structures produced by fillMissingFlows.m to calculate the average streamflow values for each HUC watershed and outputs them to the csv. Both an arithmetic average and an inverse distance weighted average are calculated for each watershed scale and a combined table is produced for each year, titled ‘flowDataFilled(YEAR).csv’



