%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The purpose of this script is to read in daily flow text files downloaded 
% from USGS (with headers removed) and store as 7-day, 14-day and 28-day 
% flows and flow percentiles for each region, in MATLAB data structures. 
% Each data structure contains an array for each site ID that is available 
% within the region. This script requires gridGaugeAlign.mat data to run. 
% The script can also calculates the percentage of missing data.

% Aug 2022
% Courtney DiVittorio
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%navigate to folder where clean text files are located (headers removed)
% list of file names to read 
finfo = dir;
fnamest = {finfo.name}';
%manually delete names that do not correspond to files
%site names - take from siteGeo structure
load('gridGaugeAlign.mat')
snames = fieldnames(siteGeo);

%%
%begin for loop for each file
for i = 1:length(fnamest)
    if i > 1
        clear percMissing30 percMissing20 percMissing7day percMissing14day percMissing28day
        tmpFile = fnamest{i};
        % removed unwanted lines in linux bash shell script so can easily read in table
        %read in table, start with arkansas
        tmptab = readtable(tmpFile,'Format','%*s %f %s %s %s %s %s %s %s %s',...
            'Delimiter', '\t', 'HeaderLines', 0,'ReadVariableNames', false);
        %if all columns filled, then want to take last streamflow value that
        %corresponds to mean (1st is max and 2nd is min)
       
        %create date vector
        dates = (datetime(1990,1,1):datetime(2022,8,1))';
        %get list of siteIDs from existing data structure
        tmpSites = siteGeo.(snames{i}).siteID;
    
        %create data structure
        flowsFullTime = struct();
        %add dates
        flowsFullTime.dates = dates;
        %for each site grab flow data and store
        k=0;
        m=0;
        tic
        for j = 1:length(tmpSites)
            %first index rows that correspond to siteID
            tmpind = find(tmptab.Var1==tmpSites(j));
            tmpdates = tmptab.Var2(tmpind);
            %make sure nothing stored in var7
            tmp = str2double(tmptab.Var7(tmpind));
            if sum(tmp)==0 || isnan(sum(tmp))==1
                tmpflow = tmptab.Var3(tmpind);
            else 
                k=k+1;
                tmpflow = tmptab.Var7(tmpind);
            end
            %need to account for missing values
            if length(tmpdates)==length(dates)
                %no missing data
               %get flow data, convert to double and store 
                tmpflow2 = str2double(tmpflow);
                %store in data structure
                flowsFullTime.(['sid',num2str(tmpSites(j))]) = tmpflow2;
            else
                %need to read in datevector and put values in correct spot
                %create empty vec
                tmpflow2 = NaN(length(dates),1);
                for n = 1:length(tmpdates)
                    tmpd = char(tmpdates{n});
                    tmpyr = str2double(tmpd(1:4));
                    tmpmo = str2double(tmpd(6:7));
                    tmpdy = str2double(tmpd(9:10));
                    dind = find(tmpyr == year(dates) & tmpmo == month(dates) & tmpdy == day(dates));
                    tmpflow2(dind) = str2double(tmpflow(n));
                end
                flowsFullTime.(['sid',num2str(tmpSites(j))]) = tmpflow2;
            end
           
        end
        toc
        disp(['part 1 region ',num2str(i),' done'])
    end
    %% calc percent of missing data over 30 yr period
    tic
    percMissing30 = NaN(length(tmpSites),1);
    %k=0;
    for j = 1:length(tmpSites)
        tmpf = flowsFullTime.(['sid',num2str(tmpSites(j))]);
    %     if nansum(tmpf) == 0   
    %         k=k+1;
    %         nansites30(k,:) = (['sid',num2str(tmpSites(j))]);
    %     end
        %calculate percentage of missing values
        percMissing30(j,1) = 100*(sum(isnan(tmpf))/length(tmpf));
    end


    %% calc percent of missing data from 2000 to present
    percMissing20 = NaN(length(tmpSites),1);
    %k=0;
    for j = 1:length(tmpSites)
        tmpf = flowsFullTime.(['sid',num2str(tmpSites(j))]);
        tmpf = tmpf(3653:end); %reduce to 2000 - 2022
    %     if nansum(tmpf) == 0   
    %         k=k+1;
    %         nansites20(k,:) = (['sid',num2str(arkSites(j))]);
    %     end
        %calculate percentage of missing values
        percMissing20(j,1) = 100*(sum(isnan(tmpf))/length(tmpf));
    end


    %% Caclulate 7-day average for each site, starting in 1990

    %want to summarize 7 days in the past on a tuesday (tuesday - monday) 
    %first tuesday to start with is Jan 9th 1990

    %create date vector for tuesdays
    firstDate = datetime(1990,1,9);
    lastDate = datetime(2022,08,01);
    dates7day=(firstDate:7:lastDate)';

    %require at least 4 data points to calculate average

    %create empty structure for 7-day flow average\
    flows7day = struct();
    flows7day.dates = dates7day; %add date vector

    for j = 1:length(tmpSites)
        tmpf = flowsFullTime.(['sid',num2str(tmpSites(j))]);
        %get 7 days of data and average until get to end of series
        flow7 = NaN(length(dates7day),1);
        for c = 1:length(dates7day)
            tmpInd = find(dates == dates7day(c));
            tmpf7 = tmpf((tmpInd-7):tmpInd-1);
            if sum(isnan(tmpf7)) < 4 %then less than 4 are NaN
                flow7(c,1) = nanmean(tmpf7);
            end
            clear tmpf7 tmpInd
        end
        %store in new structure
        flows7day.(['sid',num2str(tmpSites(j))])=flow7;
    end

    %% Calculate 14-day average for each site
    %require 8 data points

    %create date vector for tuesdays
    firstDate = datetime(1990,1,16);
    lastDate = datetime(2022,08,01);
    dates14day=(firstDate:7:lastDate)';

    %create empty structure for 7-day flow average\
    flows14day = struct();
    flows14day.dates = dates14day; %add date vector

    for j = 1:length(tmpSites)
        tmpf = flowsFullTime.(['sid',num2str(tmpSites(j))]);
        %get 7 days of data and average until get to end of series
        flow14 = NaN(length(dates14day),1);
        for c = 1:length(dates14day)
            tmpInd = find(dates == dates14day(c));
            tmpf14 = tmpf((tmpInd-14):tmpInd-1);
            if sum(isnan(tmpf14)) < 7 %then less than 7 (half) are NaN
                flow14(c,1) = nanmean(tmpf14);
            end
            clear tmpf14 tmpInd
        end
        %store in new structure
        flows14day.(['sid',num2str(tmpSites(j))])=flow14;
    end



    %% Calculate 28-day average for each site
    %require 15 datapoints
    %create date vector for tuesdays
    firstDate = datetime(1990,1,30);
    lastDate = datetime(2022,08,01);
    dates28day=(firstDate:7:lastDate)';


    %create empty structure for 7-day flow average\
    flows28day = struct();
    flows28day.dates = dates28day; %add date vector

    for j = 1:length(tmpSites)
        tmpf = flowsFullTime.(['sid',num2str(tmpSites(j))]);
        %create empty vector of flow data 
        flow28 = NaN(length(dates28day),1);
        %get 7 days of data and average until get to end of series
        for c = 1:length(dates28day)
            tmpInd = find(dates == dates28day(c));
            tmpf28 = tmpf((tmpInd-28):tmpInd-1);
            if sum(isnan(tmpf28)) < 14 %then less than 14 (half) are NaN
                flow28(c,1) = nanmean(tmpf28);
            end
            clear tmpf28 tmpInd
        end
        %store in new structure
        flows28day.(['sid',num2str(tmpSites(j))])=flow28;
    end


    %% build frequency curve for each site and each avg and covert streamflow values to percentile

    % 7 day
    %create empty structure 
    flows7dayPerc = struct();
    flows7dayPerc.dates = dates7day; %add date vector
    percMissing7day = NaN(length(tmpSites),1);
    for j = 1:length(tmpSites)
        tmpf7 = flows7day.(['sid',num2str(tmpSites(j))]);
        %need to calc number of non-nan vales for ration
        numObs = length(tmpf7) - sum(isnan(tmpf7));
        %create vector for percentile
        flow7perc = NaN(length(dates7day),1);
        tmpf7sort = sort(tmpf7); %ascending (increasing) order
        %get rid of nan values at end
        tmpf7sort(numObs+1:end) = [];
        for k = 1:length(tmpf7)
            %find row number of sorted flows that is closest to flow value to get rank
            %then calculate percent of normal
            if isnan(tmpf7(k)) == 0
                tmpInd = find(tmpf7(k) == tmpf7sort);
                % this also works [M,tmpInd] = min(abs(tmpf7sort - tmpf7(k)));
                %if more than one value, just take lowest, unless more than 2
                %then take average position
                if length(tmpInd) < 3
                    flow7perc(k,1) = (tmpInd(1)/numObs)*100;
                else 
                    %center position
                    tmpLoc = round(length(tmpInd)/2);
                    flow7perc(k,1) = (tmpInd(tmpLoc)/numObs)*100;
                end
            end
            clear tmpInd
        end
        %store in new structure
        flows7dayPerc.(['sid',num2str(tmpSites(j))])=flow7perc;
        %record %missing so can query later
        percMissing7day(j,1) = 100*(sum(isnan(flow7perc))/length(flow7perc));
    end

    %% 14 day
    %create empty structure 
    flows14dayPerc = struct();
    flows14dayPerc.dates = dates14day; %add date vector
    percMissing14day = NaN(length(tmpSites),1);
    for j = 1:length(tmpSites)
        tmpf14 = flows14day.(['sid',num2str(tmpSites(j))]);
        %need to calc number of non-nan vales for ration
        numObs = length(tmpf14) - sum(isnan(tmpf14));
        %create vector for percentile
        flow14perc = NaN(length(dates14day),1);
        tmpf14sort = sort(tmpf14); %ascending (increasing) order
        %get rid of nan values at end
        tmpf14sort(numObs+1:end) = [];
        for k = 1:length(tmpf14)
            %find row number of sorted flows that is closest to flow value to get rank
            %then calculate percent of normal
            if isnan(tmpf14(k)) == 0
                tmpInd = find(tmpf14(k) == tmpf14sort);
                % this also works [M,tmpInd] = min(abs(tmpf7sort - tmpf7(k)));
                %if more than one value, just take lowest, unless more than 2
                %then take average position
                if length(tmpInd) < 3
                    flow14perc(k,1) = (tmpInd(1)/numObs)*100;
                else 
                    %center position
                    tmpLoc = round(length(tmpInd)/2);
                    flow14perc(k,1) = (tmpInd(tmpLoc)/numObs)*100;
                end
            end
            clear tmpInd
        end
        %store in new structure
        flows14dayPerc.(['sid',num2str(tmpSites(j))])=flow14perc;
        %record %missing so can query later
        percMissing14day(j,1) = 100*(sum(isnan(flow14perc))/length(flow14perc));
    end

    %% 28 day
    %create empty structure 
    flows28dayPerc = struct();
    flows28dayPerc.dates = dates28day; %add date vector
    percMissing28day = NaN(length(tmpSites),1);
    for j = 1:length(tmpSites)
        tmpf28 = flows28day.(['sid',num2str(tmpSites(j))]);
        %need to calc number of non-nan vales for ration
        numObs = length(tmpf28) - sum(isnan(tmpf28));
        %create vector for percentile
        flow28perc = NaN(length(dates28day),1);
        tmpf28sort = sort(tmpf28); %ascending (increasing) order
        %get rid of nan values at end
        tmpf28sort(numObs+1:end) = [];
        for k = 1:length(tmpf28)
            %find row number of sorted flows that is closest to flow value to get rank
            %then calculate percent of normal
            if isnan(tmpf28(k)) == 0
                tmpInd = find(tmpf28(k) == tmpf28sort);
                % this also works [M,tmpInd] = min(abs(tmpf7sort - tmpf7(k)));
                %if more than one value, just take lowest, unless more than 2
                %then take average position
                if length(tmpInd) < 3
                    flow28perc(k,1) = (tmpInd(1)/numObs)*100;
                else 
                    %center position
                    tmpLoc = round(length(tmpInd)/2);
                    flow28perc(k,1) = (tmpInd(tmpLoc)/numObs)*100;
                end
            end
            clear tmpInd
        end
        %store in new structure
        flows28dayPerc.(['sid',num2str(tmpSites(j))])=flow28perc;
        %record %missing so can query later
        percMissing28day(j,1) = 100*(sum(isnan(flow28perc))/length(flow28perc));
    end
    toc
    disp(['part 2 region ',num2str(i),' done']);
 
    %% create names for data structures and files and save

    newFileName = [snames{i},'FlowData.mat'];

    save(newFileName,'flowsFullTime','flows7day','flows14day','flows28day',...
        'flows7dayPerc','flows14dayPerc','flows28dayPerc','gridSitesRev')


end












