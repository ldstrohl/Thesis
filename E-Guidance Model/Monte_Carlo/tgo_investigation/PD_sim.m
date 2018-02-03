% E-Guidance
% Simulation script
% Lloyd Strohl
% 07/30/17
function [results,traj] = PD_sim(IC,FC,rocket_disp,rocket_nom,options,...
    traj_record_flag,skip,seed)

rng(seed);

% % Debug
% rngtest = randi([0 10],1);
% fprintf('Inside Test: %d\n',rngtest)

%% initial conditions and rocket parameters
r0 = IC(:,1);
V0 = IC(:,2);
% dr = IC(:,3);

% % Debug
% fprintf('Inside dr: %.1f\n',dr)

% r0_nom = IC_nom(:,1);
% V0_nom = IC_nom(:,2);

rf = FC(:,1);
Vf = FC(:,2);
af = FC(:,3);

m0 = rocket_disp(1);
v_ex = rocket_disp(2);
T_max = rocket_disp(3);
T_min = rocket_disp(4);
S = rocket_disp(5);

m0_nom = rocket_nom(1);
m0_dry = rocket_nom(2);
% v_ex_nom = rocket_nom(2);
T_max_nom = rocket_nom(4);
% T_min_nom = rocket_nom(4);

guidance_law = options(1);
tgo_method = options(2);



% Mars
R_M = 3.39619*10^6;
mu = 4.282828185603917*10^13;

% Euler Angles
[pitch,yaw,roll,~] = EulerAngles(-V0,r0,V0);

%% Simulation Settings
dt = 10^-3; % time step size
dt0 = dt;
guidance_refresh_rate = 5; % Hz
tgo_tol = 0.01; % tgo update iteration tolerance
tgo_stop = 0.5; % thrust update threshold

% Compute initial tgo
g = -mu*r0/((r0'*r0)^(1.5));
switch tgo_method
    case 1
        % gravity turn
        gamma = pi/2 - acos(r0'*V0/(norm(r0)*norm(V0)));
        gm = norm(g);
        
        a_tgo = 1/gm^2;
        b_tgo = sin(gamma)*norm(V0)^2/(2*(norm(r0)-R_M)*gm^2);
        c_tgo = -(norm(V0)^2*(1+sin(gamma)^2)/(4*(norm(r0)-R_M)*gm) + 1);
        
        a_GT = (-b_tgo + sqrt(b_tgo^2 - 4*a_tgo*c_tgo))/(2*a_tgo);
        if a_GT < 0
            a_GT = (-b_tgo - sqrt(b_tgo^2 - 4*a_tgo*c_tgo))/(2*a_tgo);
        end
        tgo0 = 1.4*norm(V0)/2 * ((1+sin(gamma))/(a_GT + gm)...
            + (1-sin(gamma))/(a_GT-gm));
    case 2
        % E-Guidance iterative
        if guidance_law == 1
            throttle_estimate = 0.43;
        else
            throttle_estimate = .8;
        end
        p = 1000;
        tgo0 = 80;
        iter = 0;
        while p > tgo_tol && iter < 100
            tgo_temp = tgo0;
            dVfV = Vf-V0;
            dV = norm(dVfV + g*tgo_temp);
            tgo0 = (m0/(throttle_estimate*T_max/v_ex))*(1-exp(-dV/v_ex));
            p = norm(tgo0 - tgo_temp);
            iter = iter + 1;
        end
    case 3
        % Souza cubic
        p = [1.5*10^7+norm(g)^2/2,...
            0,...
            -2*(V0'*V0),...
            -1,...
            -12*V0'*r0 - 18*(r0'*r0)];
        tgo_roots = roots(p);
        x = ones(1,length(tgo_roots));
        for j = 1:length(tgo_roots)
            x(j) = isreal(tgo_roots(j));
        end
        p = tgo_roots(x>0);
        tgo0 = p(p>0);
    case 4
        % Apollo cubic
        % set target jerk and snap to 0
        % convert all variables to "downrange" coordinates
        downrange = rf-r0;
        downrange = [0;downrange(2);downrange(3)];
        downrange = downrange/norm(downrange);
        rf_down = rf'*downrange;
        r0_down = r0'*downrange;
        Vf_down = Vf'*downrange;
        V0_down = V0'*downrange;
        af_down = af'*downrange;
        m_term = m0-12000;
        F_term = af*m_term;
        m_dot_term = F_term/v_ex;
        p = [af_down/2,...
            (Vf_down),...
            (rf_down - r0_down)];
            
        tgo_roots = roots(p);
        x = ones(1,length(tgo_roots));
        for j = 1:length(tgo_roots)
            x(j) = isreal(tgo_roots(j));
        end
        p = tgo_roots(x>0);
        tgo0 = -p(p<0);
end
tgo0 = 67;
% Debug
fprintf('altitude: %.1f speed: %.1f tgo0: %.0f\n',norm(r0-rf),norm(V0),tgo0)

           
% navigation uncertainty
r_sig = options(3);
V_sig = options(4); % sigmas for zero mean noise

%% initialization
fig = figure('Name','Throttle vs. tgo');
fig = fig.Number;
plot_rate = options(6);
set(gca, 'XDir','reverse')
xlabel('t_{go}')
ylabel('Throttle setting')
title('Throttle vs. t_{go}')

r = r0*ones(1,(floor(tgo0)+1)/dt);
r_nav = r0;
V = V0*ones(1,(floor(tgo0)+1)/dt);
V_nav = V0;
T = ones(1,(floor(tgo0)+1)/dt);
speed = T*norm(V(:,1));
altitude = T*norm(r(:,1)) - R_M;
range = T*norm([r(2,1);r(3,1)]);
pitch = pitch*T;
yaw = yaw*T;
roll = roll*T;
alpha = T;

Lk = alpha;
Dk = alpha;
Mach = alpha;

a = zeros(3,1);
aT = a;
m_tot = zeros(1,(floor(tgo0)+1)/dt);
tgo = tgo0*ones(1,(floor(tgo0)+1)/dt);

k = 2;
refresh_counter = (1/guidance_refresh_rate)/dt;
plot_counter = floor(tgo0/dt);

%% Simulation loop
while altitude(k) > 0 && dt ~= tgo(k-1)
    
    % compute gravity vector
    g = -mu*r(:,k)/((r(:,k)'*r(:,k))^(1.5));
    
    % get navigation readings
    [r_nav,V_nav] = nav_module(r(:,k),V(:,k),r_nav,V_nav,...
        r_sig,V_sig);
    
    % update guidance periodically
    if (tgo(k) >= tgo_stop) && (refresh_counter >= (1/guidance_refresh_rate)/dt)
        
        % update commanded thrust
        throttle = guidance_module(r_nav,V_nav,g,rf,Vf,tgo(k),...
            m0_nom-m_tot(k),T_max_nom,guidance_law,af);
        T_unlimited = norm(throttle)*T_max;
        aT = throttle*T_max/(m0-m_tot(k));
        
        % thrust limiter      
        if T_unlimited > T_max
            T(:,k) = T_max;
        elseif T_unlimited < T_min
            T(:,k) = T_min;
        else
            T(:,k) = T_unlimited;
        end
        limiter = T(:,k) / T_unlimited;
        
        aT = aT*limiter;
        refresh_counter = 0;
    else
        T_unlimited = norm(aT)*(m0-m_tot(k));
        if T_unlimited > T_max
            T(:,k) = T_max;
        elseif T_unlimited < T_min
            T(:,k) = T_min;
        else
            T(:,k) = T_unlimited;
        end
    end
    
    % produce throttle plot
    if plot_counter >= floor((tgo0/dt)/plot_rate)
        throttle_plotter(r_nav,V_nav,g,rf,Vf,tgo(k),...
            m0_nom-m_tot(k),T_max_nom,guidance_law,af,fig);
        plot_counter = 0;
    else
        plot_counter = plot_counter + 1;
    end
    
    % compute mass
    m_dot = T(:,k) / v_ex;
    m_tot(k+1) = m_tot(k) + m_dot*dt;
    if m_tot(k+1) > m0-m0_dry
        T_max = 0;
        T_min = 0;
    end
    
    % inputs 
    %   tgo_method
    %   tgo(k) x
    %   tgo_tol
    %   Vf x
    %   V(:,k) x
    %   g x
    %   m0 x
    %   throttle_estimate
    %   T_max x
    %   v_ex
    %   dt
    %   k
    % outputs
    %   tgo(k+1)
    
    % update tgo
    tgo(k+1) = tgo0-dt*k;
    if tgo(k+1) < dt % force sim to end at tgo=0
        dt = tgo(k+1);
    end
    
%     if tgo_method == 2
%         p = 1000;
%         tgo(k+1) = tgo(k);
%         iter = 0;
%         while p > tgo_tol && iter < 100
%             tgo_temp = tgo(k+1);
%             dVfV = Vf-V(:,k);
%             dV = norm(dVfV + g*tgo_temp);
%             tgo(k+1) = (m0/(throttle_estimate*T_max/v_ex))*(1-exp(-dV/v_ex));
%             p = norm(tgo(k+1) - tgo_temp);
%             iter = iter + 1;
%         end
%     else
%         tgo(k+1) = tgo0-dt*k;
%     end  
%     if tgo(k+1) < dt
%         dt = tgo(k+1);
%     end
    


    % get orientation
    [pitch(k),yaw(k),roll(k),alpha(k), body_P] = EulerAngles(aT,r(:,k),V(:,k));
    
    % atmospheric model
    [L,D,Mach(k)] = atmosphere_model(altitude(k),V(:,k),alpha(k),S,body_P);
%    L = 0;
%    D = 0;
    Lk(k) = norm(L);
    Dk(k) = norm(D);
    
    % update state
    a_LD = (L+D)/ (m0 - m_tot(k));
    [r(:,k+1),V(:,k+1)] = integrator(r(:,k),V(:,k),aT+a_LD,mu,dt);
    
    speed(k+1) = norm(V(:,k+1));
    altitude(k+1) = norm(r(:,k+1)) - R_M;
    range(k+1) = norm([r(2,k+1);r(3,k+1)]);
%     range(k+1)
    
    refresh_counter = refresh_counter+1;
    k = k + 1;
end

% snip arrays to proper length
r = r(:,1:k-1);
% V = V(:,1:k-1);
T = T(1:k-1);
tgo = tgo(1:k-1);
m_tot = m_tot(1:k-1);
speed = speed(1:k-1);
altitude = altitude(1:k-1);
range = range(1:k-1);
pitch = pitch(1:k-1);
yaw = yaw(1:k-1);
roll = roll(1:k-1);
alpha = alpha(1:k-1);
Lk = Lk(1:k-1);
Dk = Dk(1:k-1);
Mach = Mach(1:k-1);
trajectory = r-rf;

%% Data recording
    if traj_record_flag
        t = dt0*[1:skip:(k-1)]';
        traj = [t,trajectory(1,1:skip:end)',trajectory(2,1:skip:end)',...
            trajectory(3,1:skip:end)',T(1:skip:end)',tgo(1:skip:end)',...
            m0-m_tot(1:skip:end)',speed(1:skip:end)',range(1:skip:end)',...
            pitch(1:skip:end)',yaw(1:skip:end)',roll(1:skip:end)',...
            alpha(1:skip:end)'];
    else
        traj = 0;
    end
    
    

%% Visualization
if options(5)
    
    
    
    figure('Name','3D Trajectory')
    plot3(trajectory(2,:),trajectory(3,:),trajectory(1,:))
    xlabel('East')
    ylabel('North')
    zlabel('Altitude')
    title('Trajectory')
    
    figure('Name','tgo estimate')
    plot(dt0*1:k-1,tgo)
    title('t_go estimate')
    xlabel('Run time (s)')
    
    figure('Name','Thrust')
    subplot(2,1,1)
    plot(dt0*1:k-1,T/T_max)
    title('Commanded Thrust')
    xlabel('Run time (s)')
    ylabel('Throttle percentage')
    subplot(2,1,2)
    plot(dt0*1:k-1,T./(9.8*(m0-m_tot)))
    title('Thrust Acceleration (Earth G)')
    ylabel('Thrust Acceleration')
    
    figure('Name','Fuel Consumption')
    plot(dt0*1:k-1,m_tot)
    title('Fuel Consumption')
    xlabel('Run time (s)')
    ylabel('Fuel Consumed')
    
    figure('Name','Speed')
    plot(dt0*1:k-1,speed)
    title('Speed')
    xlabel('Run time (s)')
    
    figure('Name','Altitude')
    plot(dt0*1:k-1,altitude)
    title('Altitude')
    xlabel('Run time (s)')
    
    figure('Name','Speed vs. Altitude')
    plot(altitude,speed)
    title('Speed vs. Altitude')
    xlabel('Altitude')
    ylabel('Speed')
    set(gca,'Xdir','reverse')
    
    figure('Name','Range vs. Altitude')
    plot(altitude,range)
    title('Range vs. Altitude')
    xlabel('Altitude')
    ylabel('Range')
    set(gca,'Xdir','reverse')
    
    figure('Name','Euler Angles')
    subplot(3,1,1)
    plot(dt0*1:k-1,pitch)
    title('Pitch')
    subplot(3,1,2)
    plot(dt0*1:k-1,yaw)
    title('Yaw')
    subplot(3,1,3)
    plot(dt0*1:k-1,roll)
    title('Roll')
    xlabel('Run time (s)')
    figure('Name','AoA')
    plot(dt0*1:k-1,alpha)
    title('Angle of Attack')
    xlabel('Run time (s)')
end

% % error checking
% figure('Name','Forces')
% subplot(3,1,1)
% title('Error checking')
% plot(dt0*1:k-1,speed)
% ylabel('Speed')
% subplot(3,1,2)
% plot(dt0*1:k-1,Lk)
% ylabel('Lift')
% subplot(3,1,3)
% plot(dt0*1:k-1,Dk)
% ylabel('Drag')
% xlabel('Run time (s)')
% 
% figure('Name','CD')
% plot(dt0*1:k-1,Mach)
% title('CD')
% xlabel('Run time (s)')
% ylabel('CD')

results = [(k-1)*dt0;m_tot(end);range(end);speed(end);...
    acosd(-aT'*g/(norm(aT)*norm(g)))];

end

