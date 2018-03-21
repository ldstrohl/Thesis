% E-Guidance
% Simulation script
% Allows rocket to coast before igniting for optimum tgo
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
threshhold = options(6);

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
tgo_stop = 0.25; % thrust update threshold

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
        tgo0 = norm(V0)/2 * ((1+sin(gamma))/(a_GT + gm)...
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

          
% navigation uncertainty
r_sig = options(3);
V_sig = options(4); % sigmas for zero mean noise

%% initialization
r = r0*ones(1,(floor(tgo0)+1)/dt);
r_nav = r0;
V = V0*ones(1,(floor(tgo0)+1)/dt);
V_nav = V0;
T = ones(1,(floor(tgo0)+1)/dt);
thrust_angle = zeros(1,(floor(tgo0)+1)/dt);
speed = T*norm(V(:,1));
altitude = T*norm(r(:,1)) - R_M;
range = T*norm([r(2,1);r(3,1)]);
pitch = pitch*T;
yaw = yaw*T;
roll = roll*T;
alpha = T;
a_GT_vis = ones(1,(floor(tgo0)+1)/dt);
sf_vis = ones(1,(floor(tgo0)+1)/dt);
aT_tracking = ones(3,(floor(tgo0)+1)/dt);

Lk = alpha;
Dk = alpha;
Mach = alpha;

a = zeros(3,1);
aT = a;
m_tot = zeros(1,(floor(tgo0)+1)/dt);
tgo = tgo0*ones(1,(floor(tgo0)+1)/dt);

refresh_counter = (1/guidance_refresh_rate)/dt;
ignition_switch = 0; % engine off for coast
range_trigger = 0;
k = 2;

%% Simulation loop
while altitude(k) > 0 && dt ~= tgo(k-1)
    
    % compute gravity vector
    g = -mu*r(:,k)/((r(:,k)'*r(:,k))^(1.5));
    
    % get navigation readings
    [r_nav,V_nav] = nav_module(r(:,k),V(:,k),r_nav,V_nav,...
        r_sig,V_sig);
    
    % update guidance periodically
    if ignition_switch
        if ((tgo(k) >= tgo_stop) && (refresh_counter >= (1/guidance_refresh_rate)/dt))
            
            % update commanded thrust
            throttle = guidance_module(r_nav,V_nav,g,rf,Vf,tgo(k),...
                m0_nom-m_tot(k),T_max_nom,guidance_law,af);
            T_unlimited = norm(throttle)*T_max;
            aT = throttle*T_max/(m0-m_tot(k));
            aT_tracking(:,k) = aT;
            
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
            aT_tracking(:,k) = aT;
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
            aT_tracking(:,k) = aT_tracking(:,k-1);
        end
        thrust_angle(k) = acosd((-V(:,k)'*throttle)/(norm(V(:,k))*norm(throttle)));
        % get orientation
        [pitch(k),yaw(k),roll(k),alpha(k), body_P] = EulerAngles(aT,r(:,k),V(:,k));
        
    else
        % check if gravity turn acceleration a_GT is greater than
        % threshhold or if covered range is greater than range to target
        % If greater, flip ignition_switch, cut wetted area S in half,
        % trigger guidance with counter, and reset tgo0
        
        % gravity turn
        gamma = pi/2 - acos(r(:,k)'*V(:,k)/(norm(r(:,k))*norm(V(:,k))));
        gm = norm(g);
        a_tgo = 1/gm^2;
        b_tgo = sin(gamma)*norm(V(:,k))^2/(2*(norm(r(:,k))-R_M)*gm^2);
        c_tgo = -(norm(V(:,k))^2*(1+sin(gamma)^2)/(4*(norm(r(:,k))-R_M)*gm) + 1);
        a_GT = (-b_tgo + sqrt(b_tgo^2 - 4*a_tgo*c_tgo))/(2*a_tgo);
        if a_GT < 0
            a_GT = (-b_tgo - sqrt(b_tgo^2 - 4*a_tgo*c_tgo))/(2*a_tgo);
        end
        V_norm = norm(V(:,k));
        sf = (V_norm^2/(2*a_GT))*cos(gamma)*(V_norm^2+2*gm*altitude(k))/(V_norm^2+gm*altitude(k))*(R_M/norm(r(:,k)));     
        
        a_GT_vis(k) = a_GT; % logging for visualization
        sf_vis(k) = sf;
        
        if a_GT >= threshhold*T_max/m0
            range_trigger = 1;
%             fprintf('acceleration\n')
        elseif sf >= threshhold*range(k)
            range_trigger = 1;
%             fprintf('range\n')
        end
%         % SPECIAL CONFIG: UNPOWERED TRAJECTORY
%         range_trigger = 0;
        
        if range_trigger
            ignition_switch = 1;
            S = S/2;
            refresh_counter = (1/guidance_refresh_rate)/dt;
            
%             % SPECIAL CONFIG: in vacuum
%             tgo(k) = 1.2*norm(V(:,k))/2 * ((1+sin(gamma))/(a_GT + gm)...
%                 + (1-sin(gamma))/(a_GT-gm));
            % SPECIAL CONFIG: in atmosphere
            tgo(k) = norm(V(:,k))/2 * ((1+sin(gamma))/(a_GT + gm)...
                + (1-sin(gamma))/(a_GT-gm));
%             
            % snip visualization logs
            a_GT_vis = a_GT_vis(1:k);
            sf_vis = sf_vis(1:k);
        end
        
        % coast
        T(:,k) = 0;
        thrust_angle(k) = 0;

        
        % orient craft for optimal angle of attack
        %   align thrust with velocity
        %   rotate down in V x r plane to achieve AoA 
        aTc = -V(:,k); 
        axis = cross(V(:,k),r(:,k));
        axis = axis/norm(axis);
        AoA = 55;
        theta = 90-AoA;
        omcd = 1-cosd(theta);
        sd = sind(theta);
        cod = cosd(theta);
        A = [cod + axis(1)^2*omcd,axis(1)*axis(2)*omcd-axis(3)*sd,axis(1)*axis(3)*omcd+axis(2)*sd;...
            axis(2)*axis(1)*omcd+axis(3)*sd,cod+axis(2)^2*omcd,axis(2)*axis(3)*omcd-axis(1)*sd;...
            axis(3)*axis(1)*omcd-axis(2)*sd,axis(3)*axis(2)*omcd+axis(1)*sd,cod+axis(3)^2*omcd];
        aTc = A'*aTc; % rotate  
        [pitch(k),yaw(k),roll(k),alpha(k), body_P] = EulerAngles(aTc,r(:,k),V(:,k));
        
    end
    
    
    % compute mass
    m_dot = T(:,k) / v_ex;
    m_tot(k+1) = m_tot(k) + m_dot*dt;
    if m_tot(k+1) > m0-m0_dry
        T_max = 0;
        T_min = 0;
    end
    
    
    % update tgo
    tgo(k+1) = tgo(k)-dt;
    if tgo(k+1) < dt % force sim to end at tgo=0
        dt = tgo(k+1);
    end
    
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
thrust_angle = thrust_angle(1:k-1);
aT_tracking = aT_tracking(:,1:k-1);

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
    plot(dt0*[1:k-1],tgo)
    title('t_go estimate')
    xlabel('Run time (s)')
    
    figure('Name','Thrust')
    subplot(2,1,1)
    plot(dt0*[1:k-1],T/T_max)
    title('Commanded Thrust')
    xlabel('Run time (s)')
    ylabel('Throttle percentage')
    subplot(2,1,2)
    plot(dt0*1:k-1,T./(9.8*(m0-m_tot)))
    title('Thrust Acceleration (Earth G)')
    ylabel('Thrust Acceleration')
    
    figure('Name','Thrust Angle')
    plot(dt0*[1:k-1],thrust_angle)
    xlabel('Run Time (s)')
    ylabel('Thrust Angle (deg)')
    title('Thrust Angle')
    
    figure('Name','Fuel Consumption')
    plot(dt0*[1:k-1],m_tot)
    title('Fuel Consumption')
    xlabel('Run time (s)')
    ylabel('Fuel Consumed')
    
    figure('Name','Speed')
    plot(dt0*[1:k-1],speed)
    title('Speed')
    xlabel('Run time (s)')
    
    figure('Name','Altitude')
    plot(dt0*[1:k-1],altitude)
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
    plot(dt0*[1:k-1],pitch)
    title('Pitch')
    subplot(3,1,2)
    plot(dt0*[1:k-1],yaw)
    title('Yaw')
    subplot(3,1,3)
    plot(dt0*[1:k-1],roll)
    title('Roll')
    xlabel('Run time (s)')
    figure('Name','AoA')
    plot(dt0*[1:k-1],alpha)
    title('Angle of Attack')
    xlabel('Run time (s)')
    
    j = dt0*[1:length(a_GT_vis)];
    figure('Name','Gravity Turn')
    yyaxis left
    plot(j,a_GT_vis,j,T_max/m0*ones(1,length(a_GT_vis)))
    ylabel('Thrust acceleration (m/s^2)')
    yyaxis right
    plot(j,sf_vis,j,range(1:length(a_GT_vis)))
%     title('Gravity Turn')
    ylabel('Downrange (m)')
    xlabel('Run time (s)')
    legend({'a_{GT}','a_{max}','s_{GT}','Downrange'},'Location','Best')
    grid on
    
    figure('Name','Thrust Vector')
    plot(dt0*[1:k-1],aT_tracking(1,:),...
        dt0*[1:k-1],aT_tracking(2,:),...
        dt0*[1:k-1],aT_tracking(3,:))
    xlabel('Run time (s)')
    ylabel('m/s^2')
    legend('Up','East','North','Location','Northwest')
end

results = [(k-1)*dt0;m_tot(end);range(end);speed(end);...
    acosd(-aT'*g/(norm(aT)*norm(g)))];

end

