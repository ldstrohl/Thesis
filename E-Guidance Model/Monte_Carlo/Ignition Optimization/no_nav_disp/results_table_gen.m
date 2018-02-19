% Read results of Monte Carlo runs and display as table
clear
close all

latex_flag = 1;
input.tableCaption = 'Comparison of Performance With and Without Navigation Error';
input.tableLabel = 'tab:navvsnonav';
input.tablePlacement = 'ht';

load('rundata_atmo.mat') % run results file

% case, run, flight_time, fuel, range, speed, angle

% initialization
scenario = unique(rundata(:,1));
flight_time = ones(length(scenario),1);
flight_time_std = flight_time;
fuel = flight_time;
fuel_std = flight_time;
range = flight_time;
range_std = flight_time;
speed = flight_time;
speed_std = flight_time;
run_count = flight_time;

% scrape and analyze
for k = 1:length(scenario)
    A = rundata(rundata(:,1)==scenario(k),:);
    flight_time(k) = mean(A(:,4));
    flight_time_std(k) = std(A(:,4));
    fuel(k) = mean(A(:,5));
    fuel_std(k) = std(A(:,5));
    range(k) = mean(A(:,6));
    range_std(k) = std(A(:,6));
    speed(k) = mean(A(:,7));
    speed_std(k) = std(A(:,7));
    run_count(k) = length(A(:,1));
end

 row1 = [run_count(end),...
        fuel(end),fuel_std(end),...
        flight_time(end),flight_time_std(end),...
        range(end),range_std(end),...
        speed(end),speed_std(end)];

load('rundata_nonav.mat') % run results file

% case, run, flight_time, fuel, range, speed, angle

% initialization
scenario = unique(rundata(:,1));
flight_time = ones(length(scenario),1);
flight_time_std = flight_time;
fuel = flight_time;
fuel_std = flight_time;
range = flight_time;
range_std = flight_time;
speed = flight_time;
speed_std = flight_time;
run_count = flight_time;

% scrape and analyze
for k = 1:length(scenario)
    A = rundata(rundata(:,1)==scenario(k),:);
    flight_time(k) = mean(A(:,4));
    flight_time_std(k) = std(A(:,4));
    fuel(k) = mean(A(:,5));
    fuel_std(k) = std(A(:,5));
    range(k) = mean(A(:,6));
    range_std(k) = std(A(:,6));
    speed(k) = mean(A(:,7));
    speed_std(k) = std(A(:,7));
    run_count(k) = length(A(:,1));
end

 row2 = [run_count(end),...
        fuel(end),fuel_std(end),...
        flight_time(end),flight_time_std(end),...
        range(end),range_std(end),...
        speed(end),speed_std(end)];



if latex_flag
    input.data = [row1;row2];
    input.tableColLabels = {'Runs','Fuel','Fuel $\sigma$','Flight Time',...
        'FT $\sigma$','Range','Range $\sigma$','Speed','Speed $\sigma$'};
    input.tableRowLabels = {'Nav Error','No Nav Error'};
    input.dataFormatMode = 'column';
    input.dataFormat = {'%1d',1,'%.1f',8};
    input.tableColumnAlignment = 'c';
    input.tableBorders = 1;
    latex = latexTable(input);
end

