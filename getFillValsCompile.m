%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The purpose of this code is to get the average streamflow values for each
% HUC watershed and output to the csv that contains the flow data with gaps
% get average and inverse distance weighted average

%Courtney Di Vittorio
%Nov 17 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%  start annual loop
for yr = 2000:2022
    %clear previous variables
    clear flowTab watershedFill7day watershedFill14day watershedFill28day
    clear avgHUC87day avgDistHUC87day avgHUC67day avgDistHUC67day avgHUC47day avgDistHUC47day
    clear avgHUC814day avgDistHUC814day avgHUC614day avgDistHUC614day avgHUC414day avgDistHUC414day
    clear avgHUC828day avgDistHUC828day avgHUC628day avgDistHUC628day avgHUC428day avgDistHUC428day
    %read in flow tables with missing values
    tmpFlowName = (['flowData',num2str(yr),'.csv']);
    flowTab = readtable(tmpFlowName);
    %read in fill values for year yr
    tmpFillData7 = ['watershedFill7day',num2str(yr),'.mat'];
    load(tmpFillData7)
    tmpFillData14 = ['watershedFill14day',num2str(yr),'.mat'];
    load(tmpFillData14)
    tmpFillData28 = ['watershedFill28day',num2str(yr),'.mat'];
    load(tmpFillData28)
    %find matching values in fill data
    %calculate average for each watershed and add to table
    %calculate inverse distance-weighted average and add to table
    %repeat for 14 and 28day
    
    %create nan arrays for each new value
    %7-day
    avgHUC87day = NaN(size(flowTab,1),1);
    avgDistHUC87day = NaN(size(flowTab,1),1);
    avgHUC67day = NaN(size(flowTab,1),1);
    avgDistHUC67day = NaN(size(flowTab,1),1);
    avgHUC47day = NaN(size(flowTab,1),1);
    avgDistHUC47day = NaN(size(flowTab,1),1);
    %14-day
    avgHUC814day = NaN(size(flowTab,1),1);
    avgDistHUC814day = NaN(size(flowTab,1),1);
    avgHUC614day = NaN(size(flowTab,1),1);
    avgDistHUC614day = NaN(size(flowTab,1),1);
    avgHUC414day = NaN(size(flowTab,1),1);
    avgDistHUC414day = NaN(size(flowTab,1),1);
    %28-day
    avgHUC828day = NaN(size(flowTab,1),1);
    avgDistHUC828day = NaN(size(flowTab,1),1);
    avgHUC628day = NaN(size(flowTab,1),1);
    avgDistHUC628day = NaN(size(flowTab,1),1);
    avgHUC428day = NaN(size(flowTab,1),1);
    avgDistHUC428day = NaN(size(flowTab,1),1);
    

    % START WITH 7 DAY
    %nanindex for 7 day flow
    nullInd = find(isnan(flowTab.percFlow7day)==1);
    tic
    for j = 1:length(nullInd)
        %Run this part if need to match dates but do not need this 
        %because j corresponds with the order of fill values considering
        %the way previous files are written
    %     %get date, grid, lon and lat
    %     tmpDate = flowTab.date(nullInd(j));
    %     %convert to date time
    %     tmpYr = round(tmpDate,-4)/10000;
    %     tmpMo = (round(tmpDate,-2) - round(tmpDate,-4))/100;
    %     tmpDay = tmpDate - round(tmpDate,-2);
    %     tmpDate2 = datetime(tmpYr,tmpMo,tmpDay);
    %     tmpGrid = flowTab.grid(nullInd(j));
    %     tmpLon = flowTab.longitude(nullInd(j));
    %     tmpLat = flowTab.latitude(nullInd(j));
    %     %find mathcing value in fill data table
    %     %where grid index aligns
    %     gridInd = find(contains(watershedFill7day.grid,char(tmpGrid)));
    %     %where date lat and lon align
    %     fillInd = find(watershedFill7day.date == tmpDate2 & ...
    %        watershedFill7day.longitude == tmpLon & ...
    %        watershedFill7day.latitude == tmpLat);
    %     % where all align
    %     fillLoc = gridInd(gridInd==fillInd);
         %huc4 7day
         tmpDist = cell2mat(watershedFill7day.huc4gaugeDist{j});
         tmpFlow = cell2mat(watershedFill7day.huc4vals{j});
         avgDistHUC47day(nullInd(j)) = invDist(tmpDist,tmpFlow);
         avgHUC47day(nullInd(j)) = mean(tmpFlow);
         %huc6 7day
         tmpDist = cell2mat(watershedFill7day.huc6gaugeDist{j});
         tmpFlow = cell2mat(watershedFill7day.huc6vals{j});
         avgDistHUC67day(nullInd(j)) = invDist(tmpDist,tmpFlow);
         avgHUC67day(nullInd(j)) = mean(tmpFlow);
         %huc8 7day
         tmpDist = cell2mat(watershedFill7day.huc8gaugeDist{j});
         tmpFlow = cell2mat(watershedFill7day.huc8vals{j});
         avgDistHUC87day(nullInd(j)) = invDist(tmpDist,tmpFlow);
         avgHUC87day(nullInd(j)) = mean(tmpFlow);
    end
    
    % 14 day
    %nanindex for 14 day flow
    nullInd = find(isnan(flowTab.percFlow14day)==1);
    for j = 1:length(nullInd)
       
         %huc4 7day
         tmpDist = cell2mat(watershedFill14day.huc4gaugeDist{j});
         tmpFlow = cell2mat(watershedFill14day.huc4vals{j});
         avgDistHUC414day(nullInd(j)) = invDist(tmpDist,tmpFlow);
         avgHUC414day(nullInd(j)) = mean(tmpFlow);
         %huc6 7day
         tmpDist = cell2mat(watershedFill14day.huc6gaugeDist{j});
         tmpFlow = cell2mat(watershedFill14day.huc6vals{j});
         avgDistHUC614day(nullInd(j)) = invDist(tmpDist,tmpFlow);
         avgHUC614day(nullInd(j)) = mean(tmpFlow);
         %huc8 7day
         tmpDist = cell2mat(watershedFill14day.huc8gaugeDist{j});
         tmpFlow = cell2mat(watershedFill14day.huc8vals{j});
         avgDistHUC814day(nullInd(j)) = invDist(tmpDist,tmpFlow);
         avgHUC814day(nullInd(j)) = mean(tmpFlow);
    end

     % 28 day
    %nanindex for 14 day flow
    nullInd = find(isnan(flowTab.percFlow28day)==1);
    for j = 1:length(nullInd)
         %huc4 7day
         tmpDist = cell2mat(watershedFill28day.huc4gaugeDist{j});
         tmpFlow = cell2mat(watershedFill28day.huc4vals{j});
         avgDistHUC428day(nullInd(j)) = invDist(tmpDist,tmpFlow);
         avgHUC428day(nullInd(j)) = mean(tmpFlow);
         %huc6 7day
         tmpDist = cell2mat(watershedFill28day.huc6gaugeDist{j});
         tmpFlow = cell2mat(watershedFill28day.huc6vals{j});
         avgDistHUC628day(nullInd(j)) = invDist(tmpDist,tmpFlow);
         avgHUC628day(nullInd(j)) = mean(tmpFlow);
         %huc8 7day
         tmpDist = cell2mat(watershedFill28day.huc8gaugeDist{j});
         tmpFlow = cell2mat(watershedFill28day.huc8vals{j});
         avgDistHUC828day(nullInd(j)) = invDist(tmpDist,tmpFlow);
         avgHUC828day(nullInd(j)) = mean(tmpFlow);
    end
    toc


    %add new variables to table and save
    flowTab = addvars(flowTab,avgHUC87day,avgDistHUC87day,...
        avgHUC67day,avgDistHUC67day,...
        avgHUC47day,avgDistHUC47day,...
        avgHUC814day,avgDistHUC814day,...
        avgHUC614day,avgDistHUC614day,...
        avgHUC414day,avgDistHUC414day,...
        avgHUC828day,avgDistHUC828day,...
        avgHUC628day,avgDistHUC628day,...
        avgHUC428day,avgDistHUC428day,'NewVariableNames',{...
        'avg.HUC8.7day','avgDist.HUC8.7day',...
        'avg.HUC6.7day','avgDist.HUC6.7day',...
        'avg.HUC4.7day','avgDist.HUC4.7day',...
        'avg.HUC8.14day','avgDist.HUC8.14day',...
        'avg.HUC6.14day','avgDist.HUC6.14day',...
        'avg.HUC4.14day','avgDist.HUC4.14day',...
        'avg.HUC8.28day','avgDist.HUC8.28day',...
        'avg.HUC6.28day','avgDist.HUC6.28day',...
        'avg.HUC4.28day','avgDist.HUC4.28day'});
    %save to csv
    newFname = (['flowDataFilled',num2str(yr),'.csv']);
    writetable(flowTab,newFname)
    
end

%% inverse distance function

function distAvg = invDist(distances,flows)
    wts = 1./distances;
    distAvg = sum(wts.*flows)/sum(wts);
end







