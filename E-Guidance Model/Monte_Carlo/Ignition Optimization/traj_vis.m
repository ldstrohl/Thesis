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
filename = 'traj_IC6.mat';
variable_name_flag = 1; % 1 to input variable names, 0 to auto generate
save_figs = 0; % 1 to generate jpgs
fignamebase = 'tgofac' % base of fig name (thesis section)

if variable_name_flag
    % input variable name strings
    var_names = {'Altitude (m)';...
        'East (m)';...
        'North (m)';...
        'Thrust (N)';...
        't_go (s)';...
        'Ship Mass (kg)';...
        'Speed (m_s)';...
        'Groundrange (m)';...
        'Pitch (deg)';...
        'Yaw (deg)';...
        'Roll (deg)';...
        'Angle of Attack (deg)'};
    %     input figure names
    fig_names = {'alt';...
        'East (m)';...
        'North (m)';...
        'thr';...
        't_go (s)';...
        'mass';...
        'spd';...
        'rng';...
        'fat';...
        'Yaw (deg)';...
        'Roll (deg)';...
        'Angle of Attack (deg)'};
else
    var_names = ''; % leave empty
end


%% File read
data = load(filename);
data = data.traj;
data = data(2:end,:);
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

% Generate figure names
fig_count = size(fig_names);
fig_names = vertcat(fig_names,strings(data_size(2)-1 - fig_count(1)-1,1));
for k = fig_count(1)+1:data_size(2)-1
    fig_names(k) = sprintf('Var%d',k);
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
        legend_entry(k) = sprintf('Case %d',k);
    end
    %     title(sprintf('%s vs. Time',var_names(j-1)))
    ylabel(sprintf('%s',var_names(j-1)))
    xlabel('Time (s)')
                legend(legend_entry)
%         legend('Atmosphere','Vacuum')
    legend('c_t = 1.0','c_t = 1.2')
%     legend('E-Guidance','APDG')
    if save_figs
        thesis_fig(gcf,strcat(sprintf('%s',fig_names(j-1)),fignamebase))
    end
    hold off
end


figure('Name','Trajectory','NumberTitle','off')
hold on
grid on
legend_entry = strings(length(traj_index)-1,1);
for k = 1:length(traj_index)-1
    E = data(traj_index(k):traj_index(k+1)-1,3);
    N = data(traj_index(k):traj_index(k+1)-1,4);
    U = data(traj_index(k):traj_index(k+1)-1,2);
    %     t = data(traj_index(k):traj_index(k+1)-1,1);
    %     v = data(traj_index(k):traj_index(k+1)-1,j);
    plot3(E,N,U)
    legend_entry(k) = sprintf('Run %d',k);
end
xlabel('East (m)')
ylabel('North (m)')
zlabel('Altitude (m)')
if save_figs
    saveas(gcf,'traj_simvsadv','fig')
end
hold off
