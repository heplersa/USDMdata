README for WEEKLY NLDAS-2 DATA

Files are stored as netCDF format 3-dimensional stuctures [LATITUDE, LONGITUDE, WEEK OF YEAR]. Each file contains data for 1 year on the modified GPCC0.5 grid 
LATITUDE: [25.25, 49.75]
LONGITUDE: [-124.75, -65.25]

WEEKS start on the first Tuesday of the year. The first week contains the data for all days within the current and previous year before the first Tuesday of the current year. The exception is in the year 2000 (see example below).
Example 1: In the year 2000, the first Tuesday occurs on Jan 4th and the first WEEK contains data for Jan 1st to Jan 3rd. The second WEEK contains data for Jan 4th to Jan 10th.
Example 2: In the year 2002, the first Tuesday occurs on Jan 1st and the first WEEK contains data for Dec 25th to Dec 31st of 2001. The second WEEK contains data for Jan 1st to Jan 7th.
Example 3: In the year 2022, the first Tuesday occurs on Jan 4th and the first WEEK contains data for Dec 28th to Jan 3rd. The second WEEK contains data for Jan 4th to Jan 10th. The final WEEK in 2022 contains data for Jun 21st to Jun 27th. This year only has 26 weeks total.

Summary of NLDAS-2 Variables:
APCP: weekly total precipitation [kg/m2], originally from Forcing File A
EVP:  weekly average evapotranspiration [kg/m2], originally from NOAH Land Surface Model
LAI: weekly average leaf area index [-], originally from NOAH Land Surface Model
PEVAP: weekly average potential evaporation [kg/m2], originally from Forcing File A
PEVPR: weekly average potential latent heat flux [W/m2], originally from NOAH Land Surface Model
SOILM: weekly average soil moisture content [kg/m2], originally from NOAH Land Surface Model
SNOD: weekly average snow depth [m], originally from NOAH Land Surface Model
SNOM: weekly average snow melt [kg/m2], originally from NOAH Land Surface Model
SNOWC: weekly average snow cover fraction [-], orginally from NOAH Land Surface Model
SSRUN: weekly average surface runoff [kg/m2], originally from NOAH Land Surface Model
TSOIL: weekly average soil temperature [K], originally from NOAH Land Surface Model
WEASD: weekly average water equivalent of accumulated snow depth [kg/m2], originally from NOAH Land Surface Model
