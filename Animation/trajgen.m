clear
close all
% Generate FlightGear 6-DoF array from single trajectory
filename = 'traj_aniv2.mat';

% % Meteor Crater
% lat_LS = 35.028257;
% lon_LS = -111.021208;
% alt_LS = 1571;

% Wadi Rum
lat_LS = 29.525665+0.1*2-0.01*2;
lon_LS = 35.444366+-.002*10-0.01*3;
alt_LS = 959-200+2;


%% File read
data = load(filename);
data = data.traj;
data = data(2:end,:);
trajlen= size(data,1);
% for k = 1:length(traj_index)-1
%     E = data(traj_index(k):traj_index(k+1)-1,3);
%     N = data(traj_index(k):traj_index(k+1)-1,4);
%     U = data(traj_index(k):traj_index(k+1)-1,2);
%     pitch = data(traj_index(k):traj_index(k+1)-1,10);
%     yaw = data(traj_index(k):traj_index(k+1)-1, 11);
%     roll = data(traj_index(k):traj_index(k+1)-1,12);
%     t = data(traj_index(k):traj_index(k+1)-1,1);
%     %     t = data(traj_index(k):traj_index(k+1)-1,1);
%     %     v = data(traj_index(k):traj_index(k+1)-1,j);
% %     plot3(E,N,U)
% %     legend_entry(k) = sprintf('Run %d',k);
% end

    E = data(1:trajlen,3);
    N = data(1:trajlen,4);
    U = data(1:trajlen,2);
    pitch = data(1:trajlen,10);
    pitch = lpf(pitch)-4.5;
    yaw = data(1:trajlen, 11);
   yaw(yaw>10*pi/180) = 10*pi/180;
   yaw = lpf(yaw);
    roll = data(1:trajlen,12);
    roll(roll<-10*pi/180) = -10*pi/180;
    roll = lpf(roll);
    t = data(1:trajlen,1);
    [lat,lon] = utm2ll(E,N,12);
    lat_del = lat-lat(end);
    lon_del = lon-lon(end);
    lat = lat_del+lat_LS;
    lon = lon_del+lon_LS;
    U = U + alt_LS;
    
    
    tdata = [t,lon,lat,U,roll,pitch,yaw];
    astfganim
    %     t = data(traj_index(k):traj_index(k+1)-1,1);
    %     v = data(traj_index(k):traj_index(k+1)-1,j);
%     plot3(E,N,U)
%     legend_entry(k) = sprintf('Run %d',k);
