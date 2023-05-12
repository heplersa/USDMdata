%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script reads in the mean flow MATLAB data structures produced by 
% getStreamflowModelGrid.m and  creates a table that can be ingested by 
% the statistical model. Each year of tabular data is exported to a csv 
% file, called flowData(YEAR).csv. A value of NaN is used when there is no 
% gauge data.

% Courtney Di Vittorio  
% 9/22/22
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'


%start by importing 7-day mean flow data
cd('G:\My Drive\Research\NSF_Drought_Prediction\streamflow_data\flow_percentages_matlab')
load('modelFlowData7day.mat', 'flows7dayPercGridMean')
%this datastructure contains a table for each model grid, that contains
%time series

%load model grid info
cd('G:\My Drive\Research\NSF_Drought_Prediction\streamflow_data')
load('gridGaugeAlign.mat')


%% make date vector 
fullDates = flows7dayPercGridMean.C6.dates;
%start jan 4th 2000
redDates = fullDates(522:end);

%for each date, get 7-day flow in grid order, and save to table 
%create string for date in yyyymmdd in 1st column, grid in 2nd, lon in 3rd,
%lat in 4th, and 7-day flow in 5th

%fieldnames for flows to align with grid order for export
flowGrids = fieldnames(flows7dayPercGridMean);
%% for each year, grab data and export to table

for yr = 2000:2022
    clear flows7day flows14day flows28day
    tic
    fnamesave = ['flowData',num2str(yr),'.csv'];
    n=0;
    clear dates lons lats grids 
    tmpdind = find(year(fullDates)==yr);
    %get vars for date, grid, lon and lat
    for j = 1:length(tmpdind)
        tmpd = fullDates(tmpdind(j));
        formatOut = 'yyyymmdd';
        tmpdstr = datestr(tmpd,formatOut);
        for k = 1:length(gridID)
            n = n+1;
            tmpg = gridID(k);
            tmplon = gridLon(k);
            tmplat = gridLat(k);
            if isempty(find(flowGrids == tmpg))== 0 %then grid is stored in data structure
                tmpflow7 = flows7dayPercGridMean.(tmpg);
                tmpflow14 = flows14dayPercGridMean.(tmpg);
                tmpflow28 = flows28dayPercGridMean.(tmpg);
                %get flow for date of interest
                tmpind7 = find(tmpflow7.dates == tmpd);
                flows7day(n,1) = tmpflow7.meanFlow(tmpind7);
                tmpind14 = find(tmpflow14.dates == tmpd);
                flows14day(n,1) = tmpflow14.meanFlow(tmpind14);
                tmpind28 = find(tmpflow28.dates == tmpd);
                flows28day(n,1) = tmpflow28.meanFlow(tmpind28);
            else
                flows7day(n,1) = NaN;
                flows14day(n,1) = NaN;
                flows28day(n,1) = NaN;
            end
            dates(n,:) = tmpdstr;
            lons(n,1) = tmplon;
            lats(n,1) = tmplat;
            grids(n,1) = tmpg;
        end
    end
    flowData = table(dates,grids,lons,lats,flows7day,flows14day,flows28day,...
    'VariableNames',{'date','grid','longitude','latitude','percFlow7day',...
        'percFlow14day','percFlow28day'});
    writetable(flowData,fnamesave);
    toc
end

