%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The purpose of this script is to average 7-day, 14-day, and 28-day
% streamflow over each model grid, and record the number of datapoints used
% in each grid over time. 2 data structures are created for the 7-day, 
% 14-day and 28-day values, one that contains all gauges within a grid and 
% one that contains the mean of all gauges withing each grid. All tables
% and arrays are stored for each model grid.

% Aug 2022
% Courtney Di Vittorio
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load data that indicates model grid for each station
load('gridGaugeAlign.mat')
load('modelGridInfo.mat')
%%
%for each model grid
%build 2D matrix of gauge data, where time is in y direction and each gauge
%is in a different column

%to find which file(s) I need to pull into MATLAB, first seartch through
%gridSites data structure and pull names of regions that contain that
%gridID. then load flow data

%then for each site ID (use gridGaugeCount to get # of sites and
%gridGaugeSites to get site ID) get vector of data

%reduce vector to 2000 - 2022 and store in table, keeping track of which
%column I am on.

%site names
snames = fieldnames(siteGeo);

%navigate to where regional flow data is stored in MATLAB data structures. (This
%data was produced using readStreamflowDataTxt.mat


%% Start with 7-day flows
%create empty data structure
flows7dayPercGrid = struct();

tic
for j= 1:length(gridID) %for each model grid
    tmpgrid = gridID(j); %model grid
    %find which regions I need
    for k = 1:length(snames) % for each region
        %get list of grids for each region
        tmpgridvec = gridSites.(snames{k});
        tmpf = find(tmpgridvec == tmpgrid);
        %if not empty, then load matlab data, grab these sites, and clear
        if isempty(tmpf)== 0
            %create empty table
            tmptab = table();
            tmpfile = [snames{k},'FlowData.mat'];
            load(tmpfile,'flows7dayPerc')
            %add dates to table
            tmptab.dates = flows7dayPerc.dates;
            %get fieldnames
            sidnames = fieldnames(flows7dayPerc);
            %get siteID's for region
            tmpsiteIDs = siteGeo.(snames{k}).siteID;
            %site IDs that fall within grid
            tmpsiteIDmatch = tmpsiteIDs(tmpf);
            for m = 1:length(tmpf) %for each siteid that falls in grid
                %create string that will match sidnames
                tmpmatch = ['sid',num2str(tmpsiteIDmatch(m,1))];
                %grab vector and store in table
                tmptab.(tmpmatch) = flows7dayPerc.(tmpmatch);
            end
            %save table under grid name
            flows7dayPercGrid.(tmpgrid) = tmptab;
        end 
        clear flows7dayPerc 
    end
    if j == 100
        toc
        tic
    end
end
toc

%% Get average 7-day flow for each grid and record number of values available for each date

%create empty data structure
flows7dayPercGridMean = struct();

%get names of grids that have streamflow data
gridsWithData = fieldnames(flows7dayPercGrid);

for j= 1:length(gridID) %for each model grid
    %grab table for grid
    tmpgrid = gridID(j); %model grid
    %need to make sure this grid has data
    if isempty(find(tmpgrid == gridsWithData)) == 0
        tmptab = flows7dayPercGrid.(tmpgrid);
        %create new table and add dates
        tmptabnew = table();
        tmptabnew.('dates') = tmptab.dates;
        tmptab = removevars(tmptab,'dates');
        %remove dates
        %convert site data to array
        tmpdata = table2array(tmptab);
        tmpnanind = isnan(tmpdata); %returns 1 if true
        %number of non nan values in each row
        tmpnumvals = size(tmpdata,2)-sum(tmpnanind,2);
        tmpMean = nanmean(tmpdata,2);
        %store date, mean, and number of values in table and save to structure
        tmptabnew.('meanFlow') = tmpMean;
        tmptabnew.('numVals') = tmpnumvals;
        %store table in structure
        flows7dayPercGridMean.(tmpgrid) = tmptabnew;
    end
end

%% Look at results

%use gridID map to place values at a point

dateLook = datetime(2011,4,19);

%create empty maps
flowMap = NaN(size(LAT,1),size(LON,1));
numValsMap = NaN(size(LAT,1),size(LON,1));
avgNumValsMap = NaN(size(LAT,1),size(LON,1));
%for each grid in map, find mean flow at this date and number of
%observations. Also calculate average number of observations

missingData = 0;
for k = 1:size(LAT,2)
    for j = 1:size(LAT,1)
        %get grid value
        tmpgrid = gridIDmap(j,k);
        %if not missing & data exists
        if ismissing(tmpgrid) == 0 && isempty(find(tmpgrid == gridsWithData)) == 0
            %get table
            tmpMeans = flows7dayPercGridMean.(tmpgrid);
            %row that matches date
            tmpind = find(tmpMeans.dates == dateLook);
            %get mean value
            tmpMean = tmpMeans.meanFlow(tmpind);
            %get number of values
            tmpVals = tmpMeans.numVals(tmpind);
            %get average number of values
            tmpValsAvg = mean(tmpMeans.numVals);
            %put in new map
            flowMap(j,k)=tmpMean;
            numValsMap(j,k)=tmpVals;
            avgNumValsMap(j,k)=tmpValsAvg;
            %calculate number of grids with nan values
            if isnan(tmpMean) == 1
                missingData = missingData+1;
            end
        end
    end
end

totMissingData = missingData+(size(gridID,1) - size(gridsWithData,1));

totPercMissing = 100*(totMissingData/size(gridID,1));

disp(['Percentage of Grids with Missing Data = ',num2str(totPercMissing)])

%% plot
close all

figure
geoshow(LAT,LON,flowMap,'DisplayType','texturemap')
colorbar
%caxis([0 20])
xlabel('Longitude')
ylabel('Latitude')
title(['7-day Percent Flows for ',datestr(dateLook)])
set(gca,'fontsize',14)
axis equal tight

figure
geoshow(LAT,LON,numValsMap,'DisplayType','texturemap')
colorbar
%caxis([0 20])
xlabel('Longitude')
ylabel('Latitude')
title(['Number of Observations in Each Grid on ',datestr(dateLook)])
set(gca,'fontsize',14)
axis equal tight
caxis([0 10])

figure
geoshow(LAT,LON,avgNumValsMap,'DisplayType','texturemap')
colorbar
%caxis([0 20])
xlabel('Longitude')
ylabel('Latitude')
title(['Average Number of Observations in Each Grid from 1990 to 2022'])
set(gca,'fontsize',14)
axis equal tight
caxis([0 10])
%% 14 day

%create empty data structure
flows14dayPercGrid = struct();
tic
for j= 1:length(gridID) %for each model grid
    tmpgrid = gridID(j); %model grid
    %find which regions I need
    for k = 1:length(snames) % for each region
        %get list of grids for each region
        tmpgridvec = gridSites.(snames{k});
        tmpf = find(tmpgridvec == tmpgrid);
        %if not empty, then load matlab data, grab these sites, and clear
        if isempty(tmpf)== 0
            %create empty table
            tmptab = table();
            tmpfile = [snames{k},'FlowData.mat'];
            load(tmpfile,'flows14dayPerc')
            %add dates to table
            tmptab.dates = flows14dayPerc.dates;
            %get fieldnames
            sidnames = fieldnames(flows14dayPerc);
            %get siteID's for region
            tmpsiteIDs = siteGeo.(snames{k}).siteID;
            %site IDs that fall within grid
            tmpsiteIDmatch = tmpsiteIDs(tmpf);
            for m = 1:length(tmpf) %for each siteid that falls in grid
                %create string that will match sidnames
                tmpmatch = ['sid',num2str(tmpsiteIDmatch(m,1))];
                %grab vector and store in table
                tmptab.(tmpmatch) = flows14dayPerc.(tmpmatch);
            end
            %save table under grid name
            flows14dayPercGrid.(tmpgrid) = tmptab;
        end 
        clear flows14dayPerc 
    end
    if j == 100
        toc
        tic
    end
end
toc


%% Get average 14-day flow for each grid and record number of values available for each date

%create empty data structure
flows14dayPercGridMean = struct();

%grids with data
%gridsWithData = fieldnames(flows14dayPercGrid);

for j= 1:length(gridID) %for each model grid
    %grab table for grid
    tmpgrid = gridID(j); %model grid
    %need to make sure this grid has data
    if isempty(find(tmpgrid == gridsWithData)) == 0
        tmptab = flows14dayPercGrid.(tmpgrid);
        %create new table and add dates
        tmptabnew = table();
        tmptabnew.('dates') = tmptab.dates;
        tmptab = removevars(tmptab,'dates');
        %remove dates
        %convert site data to array
        tmpdata = table2array(tmptab);
        tmpnanind = isnan(tmpdata); %returns 1 if true
        %number of non nan values in each row
        tmpnumvals = size(tmpdata,2)-sum(tmpnanind,2);
        tmpMean = nanmean(tmpdata,2);
        %store date, mean, and number of values in table and save to structure
        tmptabnew.('meanFlow') = tmpMean;
        tmptabnew.('numVals') = tmpnumvals;
        %store table in structure
        flows14dayPercGridMean.(tmpgrid) = tmptabnew;
    end
end

%% 28 day

%create empty data structure
flows28dayPercGrid = struct();
tic
for j= 1:length(gridID) %for each model grid
    tmpgrid = gridID(j); %model grid
    %find which regions I need
    for k = 1:length(snames) % for each region
        %get list of grids for each region
        tmpgridvec = gridSites.(snames{k});
        tmpf = find(tmpgridvec == tmpgrid);
        %if not empty, then load matlab data, grab these sites, and clear
        if isempty(tmpf)== 0
            %create empty table
            tmptab = table();
            tmpfile = [snames{k},'FlowData.mat'];
            load(tmpfile,'flows28dayPerc')
            %add dates to table
            tmptab.dates = flows28dayPerc.dates;
            %get fieldnames
            sidnames = fieldnames(flows28dayPerc);
            %get siteID's for region
            tmpsiteIDs = siteGeo.(snames{k}).siteID;
            %site IDs that fall within grid
            tmpsiteIDmatch = tmpsiteIDs(tmpf);
            for m = 1:length(tmpf) %for each siteid that falls in grid
                %create string that will match sidnames
                tmpmatch = ['sid',num2str(tmpsiteIDmatch(m,1))];
                %grab vector and store in table
                tmptab.(tmpmatch) = flows28dayPerc.(tmpmatch);
            end
            %save table under grid name
            flows28dayPercGrid.(tmpgrid) = tmptab;
        end 
        clear flows28dayPerc 
    end
    if j == 100
        toc
        tic
    end
end
toc


%% Get average 28-day flow for each grid and record number of values available for each date

%create empty data structure
flows28dayPercGridMean = struct();

%grids with data
%gridsWithData = fieldnames(flows14dayPercGrid);

for j= 1:length(gridID) %for each model grid
    %grab table for grid
    tmpgrid = gridID(j); %model grid
    %need to make sure this grid has data
    if isempty(find(tmpgrid == gridsWithData)) == 0
        tmptab = flows28dayPercGrid.(tmpgrid);
        %create new table and add dates
        tmptabnew = table();
        tmptabnew.('dates') = tmptab.dates;
        tmptab = removevars(tmptab,'dates');
        %remove dates
        %convert site data to array
        tmpdata = table2array(tmptab);
        tmpnanind = isnan(tmpdata); %returns 1 if true
        %number of non nan values in each row
        tmpnumvals = size(tmpdata,2)-sum(tmpnanind,2);
        tmpMean = nanmean(tmpdata,2);
        %store date, mean, and number of values in table and save to structure
        tmptabnew.('meanFlow') = tmpMean;
        tmptabnew.('numVals') = tmpnumvals;
        %store table in structure
        flows28dayPercGridMean.(tmpgrid) = tmptabnew;
    end
end








