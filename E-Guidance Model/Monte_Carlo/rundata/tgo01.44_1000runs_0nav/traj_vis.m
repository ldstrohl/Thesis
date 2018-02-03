% Trajectory Visualization Tool
% Plot multiple trajectories the data of which are written in a single ascii file
% continuously
% The first column must be time or another independent variable that
% increases monotically along any trajectory
%   Read provided data file
%   Plot each variable separately
%   Plot runs together
% Lloyd Strohl - 09/08/17
clear
close all

%% Inputs
filename = 'traj4.mat';
variable_name_flag = 1; % 1 to input variable names, 0 to auto generate

if variable_name_flag
    % input variable name strings
    var_names = {'Altitude (m)';...
        'East (m)';...
        'North (m)';...
        'Thrust (N)';...
        't_go (s)';...
        'Ship Mass (kg)';...
        'Speed (m/s)';...
        'Groundrange (m)';...
        'Pitch (deg)';...     
        'Yaw (deg)';...   
         'Roll (deg)';...   
         'Angle of Attack (deg)'};
else
    var_names = ''; % leave empty
end


%% File read
data = load(filename);
data = data.traj;
data_size = size(data); % (timesteps, variables)

% check variable count
if variable_name_flag
    if data_size(2)-1 < length(var_names)
        fprintf('Variables name count mismatch - data has more variables\n')
        
    elseif data_size(2)-1 > length(var_names)
        fprintf('Variable name count mismatch - data has fewer variables\n')
    end
else 
    var_names = 'Var1';
end
% Generate variable names
var_count = size(var_names);
var_names = vertcat(var_names,strings(data_size(2)-1 - var_count(1)-1,1));
for k = var_count(1)+1:data_size(2)-1
    var_names(k) = sprintf('Var%d',k);
end

% break data into run lengths
traj_index(1) = 1;
for k = 2:data_size(1)
    if data(k,1) < data(k-1,1)
        traj_index = [traj_index;k]; % traj start row
    end
end
traj_index(end+1) = data_size(1)+1;

%% Visualization
for j = 2:data_size(2)
    figure('Name',sprintf('%s',var_names(j-1)),'NumberTitle','off')
    hold on
    grid on
    legend_entry = strings(length(traj_index)-1,1);
    for k = 1:length(traj_index)-1
        t = data(traj_index(k):traj_index(k+1)-1,1);
        v = data(traj_index(k):traj_index(k+1)-1,j);
        plot(t,v)
        legend_entry(k) = sprintf('Run %d',k);
    end
    title(sprintf('%s vs. Time',var_names(j-1)))
    ylabel(sprintf('%s',var_names(j-1)))
    xlabel('Time (s)')
%     legend(legend_entry)
    hold off
end

