%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code sorts through missing flow values for each model grid and each date
% and 1) calculates the number of gauges within the HUC8,6, and 4
% watersheds that could be used to gap-fill. It also calculates the
% distance between the centroid of the model grid and the gauge locations,
% in decimal degrees and stores all data in a MATLAB structure, for each
% year

%Courtney Di Vittorio
% 9/22/22
% Updated 9/29/22 to iterate through all years and store data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%  Load data
load('gridGaugeAlign.mat')
load('modelGridInfo.mat')
load('stationHUCs.mat')
%navigate to directory where flow data that aligns with model grids is
%stored (data produced from getStreamflowModelGrid.m)
load('modelFlowData7day.mat', 'flows7dayPercGrid')
%note - change to 14day or 28day if want to run those
gaugeGeo = readtable('gaugeLocations.csv');

%% start reading through data and creating table, go year-by-year
%navigate to directory where
for yr = 2000:2022
    tic %start timer
    %load a year of flow data
    tmpflowfile = ['flowData',num2str(yr),'.csv'];
    origflows = readtable(tmpflowfile);
       
    
    %7-day flows shown here, but can modify names for 14-day and 28-day
    
    %if data is missing, get lat and lon for missing data point
    %find HUC code for lat and lon coord using gridHUCs
    %for grid to find fill values, identify 20 x 20 window of grids to search
    %using LAT & LON
    %identify stream gauges within HUC4,6,8 using stationHUCS
    %for each grid within 20 x 20 window, extract grid data structure from
    %flows7daypercgrid and search for sites that match sites within watersheds.
    %if value is present, extract values for date of interest and record site
    %id in vectors
    %for each siteID, get lat and lon coordinate and calculate distance from
    %grid centroid using gaugeGeo store distance in vector
    %for each entry, record the date, grid, and # guages with data within
    %the watershed, vector of distances to gauge, and vector of flow values.
    %repeat for each HUC
    
      
    m=0; %counter
    %if data is missing, get lat and lon for missing data point
    for d = 1:size(origflows,1)
        %find out if flow value is missing
        tmporigflow = origflows.percFlow7day(d);
        if isnan(tmporigflow)
            m = m+1;
            %get date, lat, lon and grid id
            datestr = num2str(origflows.date(d));
            datet = datetime(str2double(datestr(1:4)),str2double(datestr(5:6)),str2double(datestr(7:8)));
            gridtest = origflows.grid(d);
            longtest = origflows.longitude(d); 
            lattest = origflows.latitude(d);
            
            %find HUC code for lat and lon coord using gridHUCs
            gridlats = gridHUCs.Lat;
            gridlons = gridHUCs.Lon;
            tmpdiff = abs(gridlats - lattest)+abs(gridlons-longtest);
            tmpind = find(tmpdiff==0);
            huc4 = gridHUCs.huc4(tmpind);
            huc6 = gridHUCs.huc6(tmpind);
            huc8 = gridHUCs.huc8(tmpind);
            
            %for grid to find fill values, identify 20 x 20 window of grids to search
            %using LAT & LON
            londiff = abs(LON(1,:)-longtest);
            latdiff = abs(LAT(:,1)-lattest);
            tmprow = find(latdiff==0);
            tmpcol = find(londiff==0);
            %check 
            %gridIDmap(tmprow,tmpcol)
            %get 20x20 window, but correct if out of bounds
            rowstart = tmprow-10;
            if rowstart <= 0
                rowstart = 1;
            end
            rowend = tmprow+10;
            if rowend > size(gridIDmap,1)
                rowend = size(gridIDmap,1);
            end
            colstart = tmpcol-10;
            if colstart <= 0
                colstart = 1;
            end
            colend = tmpcol+10;
            if colend > size(gridIDmap,2)
                colend = size(gridIDmap,2);
            end
            %get gridIDs in window
            gridIDwind = gridIDmap(rowstart:rowend,colstart:colend);
            gridIDwind = gridIDwind(:);
            %get rid of missing values
            gridIDwind(ismissing(gridIDwind)==1)=[];
            
            %identify stream gauges within HUC4,6,8 that aligns with model grid using stationHUCS
            tmphuc = stationHUCs.huc4;
            tmphucind = find(tmphuc == huc4); %checked in arcmap
            %get stations for each huc
            tmpsites = stationHUCs.siteID;
            huc4sites = tmpsites(tmphucind);
            tmphuc = stationHUCs.huc6;
            tmphucind = find(tmphuc == huc6);
            huc6sites = tmpsites(tmphucind);
            tmphuc = stationHUCs.huc8;
            tmphucind = find(tmphuc == huc8);
            huc8sites = tmpsites(tmphucind);
            %no difference between 4 and 6 here
            
            %for each grid within 20 x 20 window, extract grid data structure from
            %flows7daypercgrid and search for sites that match sites within watersheds.
            %counter to add sites
            n4 = 0;
            n6 = 0;
            n8 = 0;
            clear tmphuc4fills tmphuc4dist tmphuc6fills tmphuc6dist tmphuc8fills tmphuc8dist
            for j = 1:length(gridIDwind)
                tmpgrid = gridIDwind(j);
                %if gridid does not match any fieldnames then skip b/c no data
                if isempty(find(tmpgrid == gridsWithData))
                    %do nothing
                else
                    tmpflows7 = flows7dayPercGrid.(tmpgrid);
                    tmpgridsites = tmpflows7.Properties.VariableNames;
                    %for 2nd to last entry, get siteID and convert to num
                    clear tmpgridsitesn
                    tmpgridsitesn = nan(length(tmpgridsites)-1,1);
                    for gs = 2:length(tmpgridsites)
                        tmpgs = tmpgridsites{gs};
                        tmpgridsitesn(gs-1,1)=str2double(tmpgs(4:end));
                    end
                    %convert string to number and compare number instead of string
                    %for each site add "sid" and see if existis in table
                    %huc4
                    %counter to store fill vals 
                    for k = 1:length(huc4sites)
                        tmphucsite = huc4sites(k);
                        tmptabid = ['sid',num2str(tmphucsite)];
                        %see if siteIDexists in table
                        %if it does, then grab value, but if not, move to next huc4 site
                        if isempty(find(tmphucsite == tmpgridsitesn)) == 0 %cell2mat(strfind(tmpgridsites,tmptabid))) == 0 
                            %if value found, make sure strings are the same length
                            %if value is present, extract values for date of interest and record site
                            %id in vectors
                            tmpflows7site = tmpflows7.(tmptabid);
                            tmpdind = find(tmpflows7.dates == datet);
                            tmpflow = tmpflows7site(tmpdind);
                            if isnan(tmpflow) == 0
                                %record value
                                n4 = n4+1;
                                tmphuc4fills(n4) = tmpflow;
                                %site id is tmphucsite, use this to get distance
                                tmpgind = find(gaugeGeo.siteID == tmphucsite);
                                tmpgaugelat = gaugeGeo.siteLat(tmpgind);
                                tmpgaugelon = gaugeGeo.siteLon(tmpgind);
                                %calc distance in decimel degrees
                                tmphuc4dist(n4) = ((longtest-tmpgaugelon)^2+(lattest-tmpgaugelat)^2)^0.5;
                            end
          
                        end
                    end
       
                    %repeat for huc6
                    for k = 1:length(huc6sites)
                        tmphucsite = huc6sites(k);
                        tmptabid = ['sid',num2str(tmphucsite)];
                        %see if exists in table
                        %if it does, then grab value, but if not, move to next huc4 site
                        if isempty(find(tmphucsite == tmpgridsitesn)) == 0 %isempty(cell2mat(strfind(tmpgridsites,tmptabid))) == 0 
                            tmpflows7site = tmpflows7.(tmptabid);
                            tmpdind = find(tmpflows7.dates == datet);
                            tmpflow = tmpflows7site(tmpdind);
                            if isnan(tmpflow) == 0 %value exists
                                %record value
                                n6 = n6+1;
                                tmphuc6fills(n6) = tmpflow;
                                %site id is tmphucsite, use this to get distance
                                tmpgind = find(gaugeGeo.siteID == tmphucsite);
                                tmpgaugelat = gaugeGeo.siteLat(tmpgind);
                                tmpgaugelon = gaugeGeo.siteLon(tmpgind);
                                %calc distance in decimel degrees
                                tmphuc6dist(n6) = ((longtest-tmpgaugelon)^2+(lattest-tmpgaugelat)^2)^0.5;
                            end
                        end
                    end
                    %repeat for huc8
                     for k = 1:length(huc8sites)
                        tmphucsite = huc8sites(k);
                        tmptabid = ['sid',num2str(tmphucsite)];
                        %see if exists in table
                        %if it does, then grab value, but if not, move to next huc4 site
                        if isempty(find(tmphucsite == tmpgridsitesn)) == 0 %isempty(cell2mat(strfind(tmpgridsites,tmptabid))) == 0 %~any(strcmp(tmptabid ,tmpgridsites))
                            %if value is present, extract values for date of interest and record site
                            %id in vectors
                            tmpflows7site = tmpflows7.(tmptabid);
                            tmpdind = find(tmpflows7.dates == datet);
                            tmpflow = tmpflows7site(tmpdind);
                            if isnan(tmpflow) == 0
                                %record value
                                n8 = n8+1;
                                tmphuc8fills(n8) = tmpflow;
                                %site id is tmphucsite, use this to get distance
                                tmpgind = find(gaugeGeo.siteID == tmphucsite);
                                tmpgaugelat = gaugeGeo.siteLat(tmpgind);
                                tmpgaugelon = gaugeGeo.siteLon(tmpgind);
                                %calc distance in decimel degrees
                                tmphuc8dist(n8) = ((longtest-tmpgaugelon)^2+(lattest-tmpgaugelat)^2)^0.5;
                            end
        
                        end
                     end
                end
            end
            %count number of fill values and add to vectors and cell arrays
            %if no fill values then make empty
            if n4 > 0
                huc4ns(m,1) = length(tmphuc4fills);
                huc4dist{m}= num2cell(tmphuc4dist);
                huc4vals{m}= num2cell(tmphuc4fills);
            else
                huc4ns(m,1) = 0;
                huc4dist{m}= {[NaN]};
                huc4vals{m}= {[NaN]};
            end
            if n6 > 0   
                huc6ns(m,1) = length(tmphuc6fills);
                huc6dist{m}= num2cell(tmphuc6dist);
                huc6vals{m}= num2cell(tmphuc6fills);
            else
                huc6ns(m,1) = 0;
                huc6dist{m}= {[NaN]};
                huc6vals{m}= {[NaN]};
            end
            if n8 > 0
                huc8ns(m,1) = length(tmphuc8fills); 
                huc8dist{m}= num2cell(tmphuc8dist);
                huc8vals{m}= num2cell(tmphuc8fills);
            else
                huc8ns(m,1) = 0;
                huc8dist{m}= {[NaN]};
                huc8vals{m}= {[NaN]};
            end
            grids(m,1) = gridtest;
            lons(m,1) = longtest;
            lats(m,1) = lattest;
            dates(m) = datet;
        end
    end
    
    
    %% after looping through all missing dates, store in table
    watershedFill7day = table(dates',grids,lons,lats,huc4ns,huc4dist',huc4vals',...
        huc6ns,huc6dist',huc6vals',huc8ns,huc8dist',huc8vals','VariableNames',...
        {'date','grid','longitude','latitude','huc4numFills','huc4gaugeDist','huc4vals',...
        'huc6numFills','huc6gaugeDist','huc6vals','huc8numFills','huc8gaugeDist','huc8vals'});
    toc %stop timer
    %save to matlab data file
    fname = ['watershedFill7day',num2str(yr),'.mat'];
    save(fname,"watershedFill7day")
end
