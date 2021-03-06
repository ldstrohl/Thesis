% Read results of Monte Carlo runs and display as table
clear
close all

load('rundata.mat') % run results file
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
for k = 1:length(unique(rundata(:,1)))
    A = rundata(rundata(:,1)==k,:);
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

% display table
results = table(run_count,...
    fuel,fuel_std,...
    flight_time,flight_time_std,...
    range,range_std,...
    speed,speed_std,...
    'RowNames',num2cell(num2str(scenario)),...
    'VariableNames',...
    {'Runs',...
    'Fuel', 'Fuel_dev',...
    'Flight_Time','FT_dev',...
    'Range','Range_dev',...
    'Speed','Speed_dev'})

