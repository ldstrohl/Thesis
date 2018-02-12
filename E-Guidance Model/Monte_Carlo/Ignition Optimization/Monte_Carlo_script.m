 % E-Guidance run script

clear
close all

%% Inputs
runs = 10; % per case
scenario = 6; % 1:6
debug_mode = 0;
debug_row = 12;

% Data recording
data_save_flag = 0; % save run data
result_file = 'rundata_ignopti_flip.mat';
traj_file = 'traj_ignopti_flip_'; % prefix - case is appended
traj_rate = 1/1000; % how many runs to record
traj_step = 200; % time steps per recorded line

% Messages (1 on, 0 off)
run_and_seed = 0;
run_results = 1;
case_results = 0;
vis_flag = 0; % plots of each run

% Simulation settings
guidance_law = 2; % 1 for simple, 2 for final commanded thrust
af_factor = 2; % final commanded thrust, in Mars g (only for law 2)
tgo_method = 1; % 1: GT, 2: fixed-point, 3: Souza, 4: apollo cubic
threshhold = 1; % factor of T_max required for gravity turn 
                  %     to ignite engine and start guidance

% models
rocket_dispersion_flag = 1; % 1 turns on rocket paramter dispersion, 0 for off
IC_dispersion_flag  = 1; % on/off initial condition dispersion
atmosphere_flag = 1; % on/off atmosphere model (off = vacuum)
nav_flag = 1; % on/off navigation errors


%% Conditions
% Mars
R_M = 3.39619*10^6;
mu = 4.282828185603917*10^13;

% Rocket parameters
m0_nom = 58000;
m0_dry = 1000;
T_max_nom = 800000;
T_min_nom = 200000;
v_ex_nom = 3531.7;

if atmosphere_flag == 1
    S = 62.21; % m^2 wetted area
else
    S = 0;
end

if nav_flag == 1
    r_sig = 1;
    V_sig = 1/3;
else
    r_sig = 0;
    V_sig = 0;
end

% Rocket Parameter dispersion
if rocket_dispersion_flag == 1
    rocket_disp = 2/100;
else
    rocket_disp = 0/100;
end
m0_disp = rocket_disp;
v_ex_disp = rocket_disp;
T_max_disp = rocket_disp;
T_min_disp = rocket_disp;

% Initial Condition dispersion
if IC_dispersion_flag == 1
    ICV_disp = 1/3;
    ICr_disp = 1/3;
else
    ICV_disp = 0;
    ICr_disp = 0;
end
dVmax = 10;
drmax = 1000;

% Initialization
if debug_mode
    data_save_flag = 0;
    traj_rate = 1/1000000;
end

run_fid = fopen(result_file);
if data_save_flag
    if run_fid == -1
        rundata = zeros(runs*length(scenario),8);
        % rundata(1,:) = {'Case','Run','Flight Time','Fuel','Range','Speed','Angle'};
        save(result_file,'rundata');
        prev_rows = 0;
%         s = rng;
%         save(seed_file,'s')
    else
        load(result_file)
%         load(seed_file)
        rundata = rundata(rundata(:,1)~=0,:);
        prev_rows = size(rundata);
        prev_rows = prev_rows(1);
        rundata = [rundata;zeros(runs*length(scenario),8)];
    end
else
    prev_rows = 0;
    traj_rate = 1/1000000;
end

if debug_mode
    load(result_file)
    scenario = rundata(debug_row,1);
    runs = 1;
    fprintf('DEBUG MODE\n')
end


%% Case loop
tic
q = 0; % case counter, needed if not running cases from 1
for j = scenario
    q  = q + 1;
    
    if data_save_flag
        case_file = strcat(traj_file,num2str(j),'.mat');
        case_fid = fopen(case_file);
        if case_fid == -1
            traj = zeros(1,13);
            save(case_file,'traj')
        else
            load(case_file)
        end
    end
    
    % Initial Condition
    switch j
        case 1
            r0=[-3392398.53237372;...
                -250959.063247613;...
                -9948.79949940360];
            V0=[106.359752854507;...
                127.905472084214;...
                536.969704060494];
        case 2
            r0=[-3392853.41196660;...
                -251521.394893207;...
                -12340.6407275791];
            V0=[100.641644007366;...
                127.405126758944;...
                550.159986736327];
        case 3
            r0=[-3393284.33533986;...
                -252082.100776225;...
                -14790.0944674171];
            V0=[ 94.8226733980752;...
                127.339377154286;...
                563.259809653096];
        case 4
            r0 = [-3394065.53680778;...
                -253202.004401189;...
                -19862.3977888315];
            V0 = [82.5755978749391;...
                127.252297560431;...
                589.632069823669];
        case 5
            r0=[ -3394734.98407792;...
                -254319.771535922;...
                -25169.1470738351];
            V0=[ 69.4859022032303;...
                126.547356413357;...
                616.543099680754];
        case 6
            r0=[  -3395285.92438538;...
                -255429.908623963;...
                -30715.5304665501];
            V0=[55.7770843260274;...
                125.451270269588;...
                644.130998398719];
    end
    
    % Final Condition
    lonf = 184.2;   % azimuth
    latf = 0;       % inclination
    rf = [R_M*cosd(lonf)*cosd(latf);...
        R_M*sind(lonf)*cosd(latf);...
        R_M*sind(latf)];
    Vf = (rf/norm(rf))*-1;
    af = af_factor*mu*rf/((rf'*rf)^(1.5));
    
    % coordinate transformation to landing site
    % [up;east;north]
    A = [cosd(lonf),sind(lonf),0;-sind(lonf),cosd(lonf),0;0,0,1];
    r0 = A*r0;
    V0 = A*V0;
    rf = A*rf;
    Vf = A*Vf;
    af = A*af;
    
    % Initialization
    flight_time = ones(1,runs);
    fuel = ones(1,runs);
    range_f = ones(1,runs);
    speed_f = ones(1,runs);
    angle_f = ones(1,runs);
    rng shuffle
    seed = randi([0 10000],[1,runs]); % create reproducible results
%     seed(1) = 7883;
    
    if case_results == 1
        fprintf('Case: %0d \n',j)
    end

    %% Run loop
    for k = 1:runs
        
        if run_and_seed == 1
            fprintf('Table Row: %d   Run: %d/%d\n',...
                (q-1)*runs+k+prev_rows,k,runs)
        end
        
        rng(seed(k),'twister')
        
%         if data_save_flag
%             load('seeds.mat')
%             s((q-1)*runs+k+prev_rows) = rng;
%             save('seeds.mat','s');
%         else
%             s((q-1)*runs+k+prev_rows) = rng;
%         end
        if debug_mode
%             k = rundata(debug_row,2);
            rng(rundata(debug_row,3))
        end
        
%         % Debugging
%         rngtest = randi([0 10],1);
%         fprintf('Outside Test: %d\n',rngtest)
        
        % determine if run should be recorded
        if mod(k,1/traj_rate) == 0
            traj_record_flag = 1;
        else
            traj_record_flag = 0;
        end
        
        m0 = m0_nom*(1+m0_disp*(1-2*rand(1)));
        v_ex = v_ex_nom*(1+v_ex_disp*(1-2*rand(1)));
        T_max = T_max_nom*(1+T_max_disp*(1-2*rand(1)));
        T_min = T_min_nom*(1+T_min_disp*(1-2*rand(1)));
        
        % disperse IC here
        % IC_disp = IC_nom * some dispesion method
        dV = ICV_disp*randn(3,1)*dVmax;
        dr = ICr_disp*randn(3,1)*drmax;
        r0_disp = r0 + dr;
        V0_disp = V0 + dV;
        

%         fprintf('dr outside: %.1f\n',dr)
        
        % create sim inputs
        IC = [r0_disp,V0_disp,dr];
        FC = [rf,Vf,af];
        rocket_disp = [m0;v_ex;T_max;T_min;S];
        rocket_nom = [m0_nom;m0_dry;v_ex_nom;T_max_nom;T_min_nom];
        options = [guidance_law;tgo_method;r_sig;V_sig;vis_flag;threshhold];
%         seed = rng;
        
        [sim_data,traj_add] = PD_sim(IC,FC,rocket_disp,rocket_nom,options,...
            traj_record_flag,traj_step,rng);
        
        if data_save_flag
            if traj_record_flag
                traj = [traj;traj_add];
                save(case_file,'traj');
            end
        end
        
        flight_time(k) = sim_data(1);
        fuel(k) = sim_data(2);
        range_f(k) = sim_data(3);
        speed_f(k) = sim_data(4);
        angle_f(k) = sim_data(5);
        if run_results == 1
            fprintf('Flight time: %.1f (s)\n',flight_time(k))
            fprintf('Range: %.1f (m)\n',range_f(k))
            fprintf('Speed: %.1f (m/s)\n',speed_f(k))
            fprintf('Fuel: %.0f (kg)\n\n',fuel(k))
        end
        
        if data_save_flag
            rundata((q-1)*runs+k+prev_rows,:) = [j,k,seed(k)...
                flight_time(k),fuel(k),range_f(k),speed_f(k),angle_f(k)];
            save(result_file,'rundata')
        end
    end
    
    % trim zeros from results data
    if data_save_flag
        rundata = rundata(rundata(:,1)~=0,:);
        save(result_file,'rundata')
    end
    
%     seeds(q,:) = seed;
    flight_time_mean(q) = mean(flight_time);
    fuel_mean(q) = mean(fuel);
    range_mean(q) = mean(range_f);
    speed_mean(q) = mean(speed_f);
    angle_mean(q) = mean(angle_f);
    
    flight_time_std(q) = std(flight_time);
    fuel_std(q) = std(fuel);
    range_std(q) = std(range_f);
    speed_std(q) = std(speed_f);
    angle_std(q) = std(angle_f);
    
    if case_results == 1
        fprintf('**** RESULTS ****\n')
        fprintf('Flight time: %.2f (s) \x3c3: %.2f\n',flight_time_mean(q),flight_time_std(q))
        fprintf('Fuel: %.2f (kg) \x3c3: %.2f\n',fuel_mean(q),fuel_std(q))
        fprintf('Range: %.2f (m) \x3c3: %.2f\n',range_mean(q),range_std(q))
        fprintf('Speed: %.2f (m/s) \x3c3: %.2f\n\n',speed_mean(q),speed_std(q))
    end
    
end
run_time = toc/(length(scenario)*runs);
fprintf('Time per run: %f\n',run_time)

%% Results
results.flight_time = [flight_time_mean',flight_time_std'];
results.fuel = [fuel_mean',fuel_std'];
results.range = [range_mean',range_std'];
results.speed = [speed_mean',speed_std'];
results.angle = [angle_mean',angle_std'];

results_table = table(results.fuel(:,1),results.fuel(:,2),...
    results.flight_time(:,1),results.flight_time(:,2),...
    results.range(:,1),results.range(:,2),...
    results.speed(:,1),results.speed(:,2),...
    results.angle(:,1),results.angle(:,2),...
    'RowNames',num2cell(num2str(scenario')),...
    'VariableNames',{'Fuel', 'Fuel_dev',...
    'Flight_Time','FT_dev',...
    'Range','Range_dev',...
    'Speed','Speed_dev',...
    'Angle','Angle_dev'});
%
% save('results.mat','results','seeds','results_table')
fclose('all');
