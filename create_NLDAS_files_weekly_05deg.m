%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Date Created: 19 July 2022
% Last Updtaed: 08 May 2023
% Author: Lauren Lowman
% 
% Interpolate NLDAS Data to 0.5 Degree Grid and create weekly averages for 
% all of CONUS. The weekly averages begin on the first Tuesday of the year.
% See README_NLDAS.txt file for more details on dataset created from this 
% script.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;

savefig = 0; % save figure = 1 | do not save figure = 0 
savedat = 0; % save output dataset = 1 | do not save output dataset = 0

%% INITIALIZE PARAMETERS

dtype = 'APCP'; % output data type from this script

atype = 'SUM'; % aggregation type (AVG or SUM) over weekly period

ftype = 'FORA'; % NOAH or FORA

% FORA Data
% 'APCP'    % Precipitation hourly total (kg/m2)
% 'PEVAP'   % Potential evaporation hourly total (kg/m2)

% NOAH Data
% 'EVP'     % Total evapotranspiration (kg/m2)
% 'PEVPR'   % Potential latent heat flux (W/m2)
% 'SOILM'   % Soil moisture content (kg/m2)
% 'SSRUN'	% Surface runoff (kg/m2)
% 'TSOIL'   % Soil temperature (K)
% 'LAI'     % Leaf area index (-)

% file type: forcing or mosaic
% ftype = 'FORA0125'; % FORA0125 or MOS0125 or NOAH_LSM

% temporal resolution
temp_res = 'WEEKLY'; % MONTHLY_CLIM; MONTHLY, HOURLY


%% SET PATHNAMES

pathDat = "ADD/PATH/TO/RAW/HOURLY/NLDAS/DATA/HERE/";

pathMask = "ADD/PATH/TO/NLDAS/CONUS/MASK/HERE/";

pathOut = "ADD/PATH/TO/STORE/OUTPUT/DATASET/HERE/";

%% SPATIAL GRID INFORMATION

latlim = [25 50];
lonlim = [-125 -65];

%% READ IN CONUS LAND MASK at 0.125 deg
% The NLDAS land/water mask is avaiable to download here: 
% 	https://ldas.gsfc.nasa.gov/nldas/specifications
%	https://ldas.gsfc.nasa.gov/sites/default/files/ldas/nldas/NLDAS_masks-veg-soil.nc4

finfo = ncinfo(fullfile(pathMask,'NLDAS_masks-veg-soil.nc4'));

latin = ncread(fullfile(pathMask,'NLDAS_masks-veg-soil.nc4'),'lat');
lonin = ncread(fullfile(pathMask,'NLDAS_masks-veg-soil.nc4'),'lon');

conus_mask = ncread(fullfile(pathMask,'NLDAS_masks-veg-soil.nc4'),'CONUS_mask');

% Create 2-d Matrices
[lon,lat] = meshgrid(lonin,latin);

% % store grid size
% grid_size = size(conus_mask);
% nlat = grid_size(2);
% nlon = grid_size(1);

% plot NLDAS mask
figure;
axesm('MapProjection','mercator','MapLatLim',latlim,'MapLonLim',lonlim,... 
    'MLineLocation',5,'PLineLocation',5,'MeridianLabel','on','ParallelLabel','on',...
    'MLabelParallel','south','MLabelLocation',-120:10:-70,'PLabelLocation',25:5:50,...
    'LabelRotation','on','GLineWidth',1,'GColor',[0.65 0.65 0.65],...
    'Frame','on','Grid','on','FontName','times','FontSize',10);

surfm(lat,lon,double(conus_mask'));

title('NLDAS CONUS LAND/SEA MASK');

tightmap;

%% READ IN GPCC MASK at 0.5 deg

% downloaded using WGET from: 
%	https://opendata.dwd.de/climate_environment/GPCC/full_data_monthly_v2020/05/

finfo = ncinfo(fullfile('/home/lowmanle/matlab_scripts/NSF_MPS_Drought','full_data_monthly_v2020_2011_2019_05.nc'));

LATin = ncread(fullfile('/home/lowmanle/matlab_scripts/NSF_MPS_Drought','full_data_monthly_v2020_2011_2019_05.nc'),'lat');
LONin = ncread(fullfile('/home/lowmanle/matlab_scripts/NSF_MPS_Drought','full_data_monthly_v2020_2011_2019_05.nc'),'lon');
precip = ncread(fullfile('/home/lowmanle/matlab_scripts/NSF_MPS_Drought','full_data_monthly_v2020_2011_2019_05.nc'),'precip');

% Create 2-d Matrices
[LON,LAT] = meshgrid(LONin,LATin);

% plot GPCC mask
figure;
axesm('MapProjection','mercator','MapLatLim',latlim,'MapLonLim',lonlim,... 
    'MLineLocation',5,'PLineLocation',5,'MeridianLabel','on','ParallelLabel','on',...
    'MLabelParallel','south','MLabelLocation',-120:10:-70,'PLabelLocation',25:5:50,...
    'LabelRotation','on','GLineWidth',1,'GColor',[0.65 0.65 0.65],...
    'Frame','on','Grid','on','FontName','times','FontSize',10);

surfm(LAT,LON,double(squeeze(precip(:,:,1)')));

title('GPCC Precip');

tightmap;

%% REFINE GRID TO REGION: -125,25,-65,50

LATsub = LAT;
LONsub = LON;

% determine rows of LAT that correspond to subregion
rows = find( LAT(:,1) <= 50 & LAT(:,1) >= 25 );
cols = find( LON(1,:) <= -65 & LON(1,:) >= -125 );

% subset LAT and LON
LATsub = LAT(rows,cols);
LONsub = LON(rows,cols);

% store grid size
grid_size = size(LATsub);
nlat = grid_size(1);
nlon = grid_size(2);


% plot subset
figure;
axesm('MapProjection','mercator','MapLatLim',latlim,'MapLonLim',lonlim,... 
    'MLineLocation',5,'PLineLocation',5,'MeridianLabel','on','ParallelLabel','on',...
    'MLabelParallel','south','MLabelLocation',-120:10:-70,'PLabelLocation',25:5:50,...
    'LabelRotation','on','GLineWidth',1,'GColor',[0.65 0.65 0.65],...
    'Frame','on','Grid','on','FontName','times','FontSize',10);

surfm(LATsub,LONsub,double(squeeze(precip(cols,rows,1)')));

title('GPCC Precip w Refined Grid');

tightmap;

% indsLAT = find(LAT > 50 | LAT < 25);
% indsLON = find(LON > -65 | LON <-125);
% 
% [p,q] = ind2sub(size(LAT),[indsLAT; indsLON]);

% LATsub(indsLAT) = [];
% LONsub(indsLAT) = [];

% LATsub(indsLON) = [];
% LONsub(indsLON) = [];

%% TEST READING IN RAW NLDAS DATA FILE

% Read in Forcing File A Variables

fname_FORA = 'NLDAS_FORA0125_H.A20000126.2100.002.grb.SUB.nc4';

finfo_FORA = ncinfo(fullfile(pathDat,fname_FORA));

APCP = ncread(fullfile(pathDat,fname_FORA),'APCP'); % Precipitation hourly total (kg/m2)
PEVAP = ncread(fullfile(pathDat,fname_FORA),'PEVAP'); % Potential evaporation hourly total (kg/m2)

% Read in NOAH Land Surface Model Variables

fname_NOAH = 'NLDAS_NOAH0125_H.A20000126.2100.002.grb.SUB.nc4';

finfo_NOAH = ncinfo(fullfile(pathDat,fname_NOAH));

LAI = ncread(fullfile(pathDatLAI,fname_NOAH_LAI),'LAI'); % Total evapotranspiration (kg/m2)
EVP = ncread(fullfile(pathDat,fname_NOAH),'EVP'); % Total evapotranspiration (kg/m2)
PEVPR = ncread(fullfile(pathDat,fname_NOAH),'PEVPR'); % Potential latent heat flux (W/m2)
SOILM = ncread(fullfile(pathDat,fname_NOAH),'SOILM'); % Soil moisture content (kg/m2)
depth_2 = ncread(fullfile(pathDat,fname_NOAH),'depth_2_bnds'); % measurement depths for soil moisture content (cm)
SSRUN = ncread(fullfile(pathDat,fname_NOAH),'SSRUN'); % Surface runoff (kg/m2)
TSOIL = ncread(fullfile(pathDat,fname_NOAH),'TSOIL'); % Soil temperature (K)
depth = ncread(fullfile(pathDat,fname_NOAH),'depth_bnds'); % measurement depths for soil moisture content (cm)

% Separate SOILM and TSOIL into individual layer (based on depth_2_bnds and depth_bnds variables in NetCDF file)
SOILM_0_200 = SOILM(:,:,1);
SOILM_0_100 = SOILM(:,:,2);
SOILM_0_10 = SOILM(:,:,3); % surface layer soil moisture
SOILM_10_40 = SOILM(:,:,4);
SOILM_40_100 = SOILM(:,:,5);
SOILM_100_200 = SOILM(:,:,6); 

TSOIL_0_10 = TSOIL(:,:,1);
TSOIL_10_40 = TSOIL(:,:,2);
TSOIL_40_100 = TSOIL(:,:,3);
TSOIL_100_200 = TSOIL(:,:,4);


% Plot data to test read-in
figure;
axesm('MapProjection','mercator','MapLatLim',latlim,'MapLonLim',lonlim,... 
    'MLineLocation',5,'PLineLocation',5,'MeridianLabel','on','ParallelLabel','on',...
    'MLabelParallel','south','MLabelLocation',-120:10:-70,'PLabelLocation',25:5:50,...
    'LabelRotation','on','GLineWidth',1,'GColor',[0.65 0.65 0.65],...
    'Frame','on','Grid','on','FontName','times','FontSize',10);

surfm(lat,lon,double(LAI'));
% surfm(lat,lon,double(squeeze(SOILM(:,:,6)')));
colorbar;
% caxis([0 0.35]); % EVP (kg/m2)
% caxis([0 0.35]); % 

title('NLDAS Data at 0.125 deg');

tightmap;



%% TEST UPSCALING NLDAS DATA TO 0.5 deg GPCC GRID

% Interpolate using griddata
% EVP_05 = griddata(lon,lat,EVP',LONsub,LATsub);
PEVPR_05 = griddata(lon,lat,PEVPR',LONsub,LATsub);
LAI_05 = griddata(lon,lat,LAI',LONsub,LATsub);


% Plot data to test interpolation
figure;
axesm('MapProjection','mercator','MapLatLim',latlim,'MapLonLim',lonlim,... 
    'MLineLocation',5,'PLineLocation',5,'MeridianLabel','on','ParallelLabel','on',...
    'MLabelParallel','south','MLabelLocation',-120:10:-70,'PLabelLocation',25:5:50,...
    'LabelRotation','on','GLineWidth',1,'GColor',[0.65 0.65 0.65],...
    'Frame','on','Grid','on','FontName','times','FontSize',10);

surfm(LATsub,LONsub,double(LAI_05));
colorbar;
% caxis([0 0.35]); % EVP (kg/m2)

title('NLDAS Data at 0.5 deg');

tightmap;


%% LOOP OVER ALL YEARS

for yyyy = 2000:2022 %2000:2022

%% DETERMINE FIRST TUESDAY OF CURRENT YEAR AND LAST TUESDAY OF PRIOR YEAR

year = yyyy;

firstTues = nweekdate(1,3,year,1); % determine date of first Tuesday of year

firstTuesDV = datevec(firstTues); % convert datenum to datevec

et = firstTuesDV(3)-1; % DOY for first Tuesday


%% CREATE LIST OF FILES FOR PREVIOUS YEAR

% grab files for days preceding first Tuesday in previous year (year 2001+)
if( year > 2000 && et < 7) % if a full week from 1st Tues starts a year prior    
    
    lastTues = lweekdate(3, year - 1, 12); % determine last Tuesday in previous year
    
    lastTuesDV = datevec(lastTues); % convert datenum to datevec
    
  
    % determine # of days needed
    nd = 7 - et;
    
    % read in last 'nd' days of previous year
    ndp = yeardays(year - 1); % # ndays in previous year

    % determine first day to read in from previous year
    st = ndp - nd;
    

    % list of all NOAH files for selected year
    files_p = dir(fullfile(pathDat,['*' ftype '*A' num2str(year-1) '*.nc4']));


    % grab only files for days preceding first Tuesday
    ind_p = st*24 + 1:length(files_p); % initial indices, updated in loop for each week
    
    % create empty matrix to store data
    dat_p = zeros(size(LATsub));

    % loop over each file to extract data
    for k = ind_p
    
        % read in file
        testNoah = ncread(fullfile(pathDat,files_p(k).name),dtype); 
        
        % Separate SOILM and TSOIL into individual layer
        if( strcmp(dtype,'SOILM') )
            testNoah = testNoah(:,:,3); % 0 - 10 cm (surface layer)
        elseif( strcmp(dtype,'TSOIL') )
            testNoah = testNoah(:,:,1); % 0 - 10 cm (surface layer)
        end
    
        % Interpolate to 0.5 deg grid using griddata
        int_dat_p = griddata(lon,lat,testNoah',LONsub,LATsub);
    
        % add to existing values (sum)
        dat_p = dat_p + int_dat_p;
       
    end
    
    disp(['First File Read in for Week 1']);
    disp(files_p(ind_p(1)).name);
    disp(['Last File Read in for Week 1']);
    disp(files_p(k).name);
 
end

%% CREATE LIST OF FILES FOR CURRENT YEAR

% list of all NOAH files for selected year
files = dir(fullfile(pathDat,['*' ftype '*A' num2str(year) '*.nc4']));

% determine # of weeks in year
if(year == 2022)
    w = 26; % until 1st week of July
else
    w = week(datetime(year,12,31));
end

% create tensor to store weekly averages
sz = size(LATsub);
out_dat = zeros(sz(1),sz(2),w); % [lat, lon, # weeks]

% grab only files for days preceding first Tuesday
inds = 1:et*24; % initial indices, updated in loop for each week


% loop over each week in the year
tic
for i = 1:w
    
    % create empty matrix to store data
    dat = zeros(size(LATsub));

    % loop over each file to extract data
    for m = inds
    
        % read in file
        testNoah = ncread(fullfile(pathDat,files(m).name),dtype);

        
        % Separate SOILM and TSOIL into individual layer
        if( strcmp(dtype,'SOILM') )
            testNoah = testNoah(:,:,3); % 0 - 10 cm (surface layer)
        elseif( strcmp(dtype,'TSOIL') )
            testNoah = testNoah(:,:,1); % 0 - 10 cm (surface layer)
        end
    
        % Interpolate to 0.5 deg grid using griddata
        int_dat = griddata(lon,lat,testNoah',LONsub,LATsub);
    
        % add to existing values (sum)
        dat = dat + int_dat;
                    
    end
        
    if( i == 1 && year > 2000 ) % if in the first week, add previous year data
        dat = dat + dat_p;
    end
    
    if( isempty(inds) )
        m = 0;
    else
        disp(['First File Read in for Week ' num2str(i)]);
        disp(files(inds(1)).name);
        disp(['Last File Read in for Week ' num2str(i)]);
        disp(files(m).name);
    end

    % calculate sum and avg of dat for past week
    if ( i == 1 && year > 2000 )
        dat_avg = dat./length([ind_p inds]);
        dat_sum = dat;
    else
        dat_avg = dat./length(inds);
        dat_sum = dat;
    end
    
    if( strcmp(atype,'AVG') )
        % save weekly averages to tensor
        out_dat(:,:,i) = dat_avg;
    else
        % save weekly accumulations to tensor
        out_dat(:,:,i) = dat_sum;
    end
    
    % update inds for next week Tues - Mon, inclusive
    if( (m+1)+7*24-1 > yeardays(year)*24 ) % if at end of year before week ends
        inds = m+1:yeardays(year)*24;
    else
        inds = m+1:(m+1)+7*24-1;
    end
    
    if( isempty(inds) && m > 0)
        break
    end

end
toc

%% Plot data to test output tensor
figure;
axesm('MapProjection','mercator','MapLatLim',latlim,'MapLonLim',lonlim,... 
    'MLineLocation',5,'PLineLocation',5,'MeridianLabel','on','ParallelLabel','on',...
    'MLabelParallel','south','MLabelLocation',-120:10:-70,'PLabelLocation',25:5:50,...
    'LabelRotation','on','GLineWidth',1,'GColor',[0.65 0.65 0.65],...
    'Frame','on','Grid','on','FontName','times','FontSize',10);

surfm(LATsub,LONsub,squeeze(out_dat(:,:,25)));
colorbar;
% caxis([0 0.35]);

title('Output Weekly Avg Tensor Data');

tightmap;

%% SAVE NETCDF FILE

if savedat

%   create a filename for dataset
    fname_out = ['NLDAS_' dtype '_' atype '_' num2str(year) '.nc'];

%     Write out data
    nccreate(fullfile(pathOut,dtype,fname_out),[dtype '_' atype],'Dimensions',{'lat',nlat,'lon',nlon,'week of year',w});
    ncwrite(fullfile(pathOut,dtype,fname_out),[dtype '_' atype],out_dat);   
    
    if( strcmp(dtype,'EVP') )
        ncwriteatt(fullfile(pathOut,dtype,fname_out),[dtype '_' atype],'Description','Weekly average total evapotranspiration');            
        ncwriteatt(fullfile(pathOut,dtype,fname_out),[dtype '_' atype],'Units','kg/m2');   
    elseif( strcmp(dtype,'PEVPR') ) 
        ncwriteatt(fullfile(pathOut,dtype,fname_out),[dtype '_' atype],'Description','Weekly average potential latent heat flux');            
        ncwriteatt(fullfile(pathOut,dtype,fname_out),[dtype '_' atype],'Units','W/m2');   
    elseif( strcmp(dtype,'SOILM') )
        ncwriteatt(fullfile(pathOut,dtype,fname_out),[dtype '_' atype],'Description','Weekly average soil moisture content');            
        ncwriteatt(fullfile(pathOut,dtype,fname_out),[dtype '_' atype],'Units','kg/m2');   
    elseif( strcmp(dtype,'SSRUN') )
        ncwriteatt(fullfile(pathOut,dtype,fname_out),[dtype '_' atype],'Description','Weekly average surface runoff');            
        ncwriteatt(fullfile(pathOut,dtype,fname_out),[dtype '_' atype],'Units','kg/m2');   
    elseif( strcmp(dtype,'TSOIL') )
        ncwriteatt(fullfile(pathOut,dtype,fname_out),[dtype '_' atype],'Description','Weekly average soil temperature');            
        ncwriteatt(fullfile(pathOut,dtype,fname_out),[dtype '_' atype],'Units','K');   
    elseif( strcmp(dtype,'APCP') )
        ncwriteatt(fullfile(pathOut,dtype,fname_out),[dtype '_' atype],'Description','Weekly average total precipitation');            
        ncwriteatt(fullfile(pathOut,dtype,fname_out),[dtype '_' atype],'Units','kg/m2');      
    elseif( strcmp(dtype,'PEVAP') )
        ncwriteatt(fullfile(pathOut,dtype,fname_out),[dtype '_' atype],'Description','Weekly average total potential evaporation');            
        ncwriteatt(fullfile(pathOut,dtype,fname_out),[dtype '_' atype],'Units','kg/m2');
    elseif( strcmp(dtype,'LAI') )
        ncwriteatt(fullfile(pathOut,dtype,fname_out),[dtype '_' atype],'Description','Weekly average leaf area index (LAI)');            
        ncwriteatt(fullfile(pathOut,dtype,fname_out),[dtype '_' atype],'Units','m2/m2');
    end

    % Write out LAT and LON
    nccreate(fullfile(pathOut,dtype,fname_out),'lat','Dimensions',{'lat',1,nlat});
    ncwrite(fullfile(pathOut,dtype,fname_out),'lat',unique(LATsub));         
    nccreate(fullfile(pathOut,dtype,fname_out),'lon','Dimensions',{'lon',1,nlon});
    ncwrite(fullfile(pathOut,dtype,fname_out),'lon',unique(LONsub));   

end

%% END LOOP
end