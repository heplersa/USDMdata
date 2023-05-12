%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Courtney Di Vittorio 
% Aug 4th 2022

% This script reads in lat and lon text files for stream gauges (downloaded 
% from USGS and finds the closest model grid and assigns a value. 
% The results from this are saved in the gridGaugeAlign.mat file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load model grid info that will be matched to gauge stations
load('modelGridInfo.mat')

%navigate to folder with lat/lon txt files

%get file names to iterate through
finfo = dir;
fnames = {finfo.name}';
%note 1st 2 rows should be skipped

%% open files, remove 1st 26 rows, and save files to new folder
for k = 3:length(fnames)
    tmpf = fnames{k};
    %record 1st 10 characters to name final vector
    fnadd = tmpf(1:10);
    %open file, remove lines, and save to file in new folder
    fid = fopen(tmpf, 'r') ;    
    for j = 1:26 %number of rows to remove
        fgetl(fid); %removes top line
    end
    buffer = fread(fid, Inf) ; % Read rest of the file.
    fclose(fid);
    %navigate to folder where you want to store cleaned text files
    fid = fopen(tmpf, 'w')  ;   % Open destination file.
    fwrite(fid, buffer) ;       % Save to file.
    fclose('all') ;

end

%% from new folder with cleaned text files, read txt into table, find
% closest model grid, and record, also save site info in table.

%need to create strings for saving lat/lon and grid info
sitenms = {'arkansas','california','greatBasin','greatLakes',...
    'lowerColorado','lowerMississippi','midAtlantic','missouri','newEngland',...
    'ohio','pacificNW','rioGrande','sourisRed','southAtlantic','tennessee',...
    'texasGulf','upperColorado','upperMississippi'};

%create data structure for each
siteGeo = struct;
gridSites = struct;
%for each file, open, remove lines, and save to file in a new folder
for k = 3:length(fnames)
   tmpf = fnames{k};
   %record 1st 10 characters to name final vector
   fnadd = tmpf(1:10); 
   tmptab = readtable(tmpf,'Format','%f %f %f %*s %*s');
   tmptab.Properties.VariableNames{'Var1'}= 'siteID';
   tmptab.Properties.VariableNames{'Var2'}= 'siteLat';
   tmptab.Properties.VariableNames{'Var3'}= 'siteLon';
   %For each station, find closes grid
   clear closestGrid
   for j = 1:size(tmptab,1)
       %grab lat and lon
       tmplat = tmptab.siteLat(j);
       tmplon = tmptab.siteLon(j);
       %tot difference between lat and lon
       geoDiff = abs(gridLat-tmplat)+abs(gridLon-tmplon);
       [minval,minloc] = min(geoDiff);
       closestGrid(j,1) = gridID(minloc);
       clear tmplat tmplon
   end
   %save with name for region added
   gridSites.(sitenms{k-2}) = closestGrid;
   siteGeo.(sitenms{k-2}) = tmptab;
   clear tmptab
end

%% For each grid go through gauge labels, record site names in columns, and count total number of gauges for each grid

%vector for total number of stream gauges
gridGaugeCount = zeros(length(gridID),1);

%maxtrix of strings that shows site numbers for each grid
gridGaugeSites =zeros(length(gridID),1);

for j = 1:length(gridID)
    tmpid = gridID(j);
    for k = 1:length(sitenms)
        %grids that match gauges
        tmpgs = gridSites.(sitenms{k});
        %siteIDs
        tmpsid = siteGeo.(sitenms{k});
        tmpsid2 = tmpsid.siteID;
        %locate where stream gauge grid aligns with model grid
        tmpmatch = find(tmpid == tmpgs);
        %if empty then move on to next region
        if isempty(tmpmatch)==0
            curcount = gridGaugeCount(j);
            %count locations and add to vector
            gridGaugeCount(j) = gridGaugeCount(j)+length(tmpmatch);
            %record site names in columns
            for m = 1:length(tmpmatch)
                gridGaugeSites(j,curcount+m) = (tmpsid2(tmpmatch(m)));
            end
        end
    end
end




%% use lat and lon to plot number of gauges for each grid

%make meshgrid for lat and lon
x = (min(gridLon):0.5:max(gridLon));
y = (max(gridLat):-0.5:min(gridLat));
[LON, LAT] = meshgrid(x,y);

gaugeCountMap = NaN(size(LON,1),size(LON,2));
%find grid ID associated with each lat and lon (if any), then find number
%of stations, and store in matrix

for k = 1:length(x)
    for j = 1:length(y)
        tmplon = LON(j,k); 
        tmplat = LAT(j,k);
        tmpind = find((gridLat == tmplat) & (gridLon == tmplon));
        if isempty(tmpind)==0
            %record gridid
            gridIDmap(j,k) = gridID(tmpind);
            %record number of gauges
            gaugeCountMap(j,k) = gridGaugeCount(tmpind);
        end
    end
end

mesh(LON,LAT,gaugeCountMap)

imagesc(gaugeCountMap)
colorbar
caxis([0 20])
axis equal tight

figure
geoshow(LAT,LON,gaugeCountMap,'DisplayType','texturemap')
colorbar
caxis([0 20])
xlabel('Longitude')
ylabel('Latitude')
title('Number of Stream Gauges for Each Model Grid')
set(gca,'fontsize',14)
sum(gaugeCountMap(:) == 0)
%448
sum(gaugeCountMap(:) > 0)
%2891
