% Read results of Monte Carlo runs and display as table
clear
close all

latex_flag = 0;
histogram_flag = 1;

fig_base = 'disppowatmo';

input.tableCaption = 'Performance of APDG In Vacuum';
input.tableLabel = 'tab:disppowatmo';
input.tablePlacement = 'ht';

load('rundata_atmo.mat') % run results file
% case, run, flight_time, fuel, range, speed, angle


% initialization
scenario = unique(rundata(:,1));
flight_time = ones(length(scenario),1);
flight_time_std = flight_time;
    flight_time_min = flight_time;
    flight_time_max = flight_time;
fuel = flight_time;
fuel_std = flight_time;
fuel_max = flight_time;
fuel_min = flight_time;
range = flight_time;
range_std = flight_time;
range_max = flight_time;
range_min = flight_time;
speed = flight_time;
speed_std = flight_time;
speed_max = flight_time;
speed_min = flight_time;
run_count = flight_time;

% scrape and analyze
for k = 1:length(scenario)
    A = rundata(rundata(:,1)==scenario(k),:);
    flight_time(k) = mean(A(:,4));
    flight_time_min(k) = min(A(:,4));
    flight_time_max(k) = max(A(:,4));
    flight_time_std(k) = std(A(:,4));
    fuel(k) = mean(A(:,5));
    fuel_min(k) = min(A(:,5));
    fuel_max(k) = max(A(:,5));
    fuel_std(k) = std(A(:,5));
    range(k) = mean(A(:,6));
    range_min(k) = min(A(:,6));
    range_max(k) = max(A(:,6));
    range_std(k) = std(A(:,6));
    speed(k) = mean(A(:,7));
    speed_min(k) = min(A(:,7));
    speed_max(k) = max(A(:,7));
    speed_std(k) = std(A(:,7));
    run_count(k) = length(A(:,1));
    if histogram_flag
        % Fuel
        h = figure('Name',strcat('Case: ',num2str(scenario(k))));
        subplot(2,1,1)
        histogram(A(:,5))
        xlabel('Fuel (kg)')
        ylabel('Run Count')
        subplot(2,1,2)
        histogram(A(:,7))
        xlabel('Speed (m/s)')
        ylabel('Run Count')
%         thesis_fig(h,strcat('hist',num2str(scenario(k)),fig_base))
    end
end

% % display table
% results = table(run_count,...
%     fuel,fuel_std,...
%     flight_time,flight_time_std,...
%     range,range_std,...
%     speed,speed_std,...
%     'RowNames',num2cell(num2str(scenario)),...
%     'VariableNames',...
%     {'Runs',...
%     'Fuel', 'Fuel_dev',...
%     'Flight_Time','FT_dev',...
%     'Range','Range_dev',...
%     'Speed','Speed_dev'})

if latex_flag
    input.data = [run_count,...
        fuel,fuel_std,fuel_max,fuel_min,...
        flight_time,flight_time_std,flight_time_max,flight_time_min,...
        range,range_std,range_max,range_min...
        speed,speed_std,speed_max,speed_min]';
    input.tableRowLabels = {'Runs','Fuel (k)','Fuel $\sigma$','Fuel max','Fuel min','Flight Time (s)',...
        'FT $\sigma$','FT max','FT min','Range (m)','Range $\sigma$','Range max','Range min','Speed (m/s)','Speed $\sigma$','Speed max','Speed min'}';
    input.tableColLabels = num2cell(num2str(scenario))';
    input.dataFormatMode = 'row';
    input.dataFormat = {'%1d',1,'%.1f',16};
    input.tableColumnAlignment = 'c';
    input.tableBorders = 1;
    latex = latexTable(input);
end
