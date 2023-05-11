# USDMdata
Discretized US Drought Data to Support Statistical Modeling

The complete data set can be found at https://github.com/heplersa/USDMdata

This repository contains R and Matlab scripts that were used to generate a discretized data set on United State drought status with other relevant environmental variables.

RasterizeUSDMforDryad.R is a R script that imports shapefiles of drought from https://droughtmonitor.unl.edu/DmData/GISData.aspx, and rasterizes them to a 0.5 degree grid and saves a flat file csv.

TeleconnectionsForDryad.R is an R script that imports raw data on four teleconnections (PNA, ENSO, NAO, AO) from cpc.ncep.noaa.gov and converts them into weekly data and exports a flat file csv.

create_NLDAS_files_weekly_05deg.m is a Matlab script that imports netCDF data files from https://disc.gsfc.nasa.gov/datasets?keywords=NLDAS and rasterizes the variables to a weekly 0.5 degree grid.
