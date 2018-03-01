% Read results of Monte Carlo runs and display as table
clear
close all

latex_flag = 0;
histogram_flag = 1;

fig_base = 'navvsnonav';

input.tableCaption = 'Performance of APDG With and Without Navigation Error';
input.tableLabel = 'tab:navvsnonav';
input.tablePlacement = 'ht';

load('rundata_nonav.mat') % run results file
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

if histogram_flag
    h_fuel = figure('Name','Fuel');
    set(gcf,'pos',[10 10 600 1500])
    h_spd = figure('Name','Speed');
    set(gcf,'pos',[10 10 600 1500])
    h_rng = figure('Name','Range');
    set(gcf,'pos',[10 10 600 1500])
end

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
    if histogram_flag && scenario(k) == 7
        % Fuel
        figure(h_fuel);
        subplot(2,1,1)
        histogram(A(:,5),[8900:100:10700])
        ylabel('Nominal')
        axis_resize(gca)
        set(gca,'xtick',[])
        
        % Range
        figure(h_rng)
        subplot(2,1,1)
        histogram(A(:,6),[0:.5:9])
        ylabel('Nominal')
        axis_resize(gca)
        set(gca,'xtick',[])
        
        % Speed
        figure(h_spd)
        subplot(2,1,1)
        histogram(A(:,7),[0:.5:15.5])
        ylabel('Nominal')
        axis_resize(gca)
        set(gca,'xtick',[])
%         thesis_fig(h,strcat('histE',fig_base))
    end
end

row1 = [run_count(end),...
        fuel(end),fuel_std(end),fuel_max(end),fuel_min(end),...
        flight_time(end),flight_time_std(end),flight_time_max(end),flight_time_min(end),...
        range(end),range_std(end),range_max(end),range_min(end)...
        speed(end),speed_std(end),speed_max(end),speed_min(end)];


load('rundata_atmo.mat')

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
if histogram_flag && scenario(k) == 7
        % Fuel
        figure(h_fuel);
        subplot(2,1,2)
        histogram(A(:,5),[8900:100:10700])
        ylabel('Dispersed')
        xlabel('Fuel (kg)')
        axis_resize(gca)
        saveas(h_fuel,strcat('hfuel',fig_base),'pdf')
%         set(gca,'xtick',[])
        
        % Range
        figure(h_rng)
        subplot(2,1,2)
        histogram(A(:,6),[0:.5:9])
        ylabel('Dispersed')
        xlabel('Range (m)')
        axis_resize(gca)
        saveas(h_rng,strcat('hrng',fig_base),'pdf')
%         set(gca,'xtick',[])
        
        % Speed
        figure(h_spd)
        subplot(2,1,2)
        histogram(A(:,7),[0:.5:15.5])
        ylabel('Dispersed')
        xlabel('Speed (m/s)')
        axis_resize(gca)
        saveas(h_spd,strcat('hspd',fig_base),'pdf')
%         set(gca,'xtick',[])
%         thesis_fig(h,strcat('histE',fig_base))
    end
end

row2 = [run_count(end),...
        fuel(end),fuel_std(end),fuel_max(end),fuel_min(end),...
        flight_time(end),flight_time_std(end),flight_time_max(end),flight_time_min(end),...
        range(end),range_std(end),range_max(end),range_min(end)...
        speed(end),speed_std(end),speed_max(end),speed_min(end)];

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
    input.data = [row1;row2]';
    input.tableRowLabels = {'Runs','Fuel (k)','Fuel $\sigma$','Fuel max','Fuel min','Flight Time (s)',...
        'FT $\sigma$','FT max','FT min','Range (m)','Range $\sigma$','Range max','Range min','Speed (m/s)','Speed $\sigma$','Speed max','Speed min'}';
    input.tableColLabels = {'Nav Error','No Nav Error'};
    input.dataFormatMode = 'row';
    input.dataFormat = {'%1d',1,'%.1f',16};
    input.tableColumnAlignment = 'c';
    input.tableBorders = 1;
    latex = latexTable(input);
end
