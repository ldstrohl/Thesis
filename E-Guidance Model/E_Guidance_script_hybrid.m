% E-Guidance
% Hybrid tgo - GT + fixed point
% Lloyd Strohl
% 5/21/17

clear
close all

%% Conditions
% Mars
R_M = 3.39619*10^6;
mu = 4.282828185603917*10^13;

% Initial Condition

for scenario = [1]
    switch scenario
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
    
    % coordinate transformation to landing site
    A = [cosd(lonf),sind(lonf),0;-sind(lonf),cosd(lonf),0;0,0,1];
    r0 = A*r0;
    V0 = A*V0;
    rf = A*rf;
    Vf = A*Vf;
    
    % Rocket parameters
    m0 = 58000;
    T_max = 800000;
    T_min = 200000;
    v_ex = 3531.7;
    
    % gravity turn tgo estimate
    gamma = pi/2 - acos(r0'*V0/(norm(r0)*norm(V0)));
    gm = norm(mu*r0/((r0'*r0)^(1.5)));
    
    a_tgo = 1/gm^2;
    b_tgo = sin(gamma)*norm(V0)^2/(2*(norm(r0)-R_M)*gm^2);
    c_tgo = -(norm(V0)^2*(1+sin(gamma)^2)/(4*(norm(r0)-R_M)*gm) + 1);
    
    a_GT = (-b_tgo + sqrt(b_tgo^2 - 4*a_tgo*c_tgo))/(2*a_tgo);
    if a_GT < 0
        a_GT = (-b_tgo - sqrt(b_tgo^2 - 4*a_tgo*c_tgo))/(2*a_tgo);
    end
    
    tgo0 = round(norm(V0)/2 * ((1+sin(gamma))/(a_GT + gm)...
        + (1-sin(gamma))/(a_GT-gm)),0);
    
    
    
    %% Simulation Settings
    dt = 10^-3; % step size
    tgo_tol = 0.01;
    % thrust_estimate = (T_max+T_min)/(T_max*2);
    thrust_estimate = norm(a_GT)*m0/T_max;
    
    
    % initialization
    r = r0*ones(1,tgo0/dt);
    V = V0*ones(1,tgo0/dt);
    T = ones(1,tgo0/dt);
    speed = T*norm(V(:,1));
    altitude = T*norm(r(:,1)) - R_M;
    range = T*norm([r(2,1);r(3,1)]);
    thrust_estimate = T*thrust_estimate;
    
    
    a = zeros(3,1);
    aT = a;
    m_tot = zeros(1,tgo0/dt);
    tgo = tgo0*ones(1,tgo0/dt);
    
    t = 0;
    k = 1;
    
    %% Simulation loop
    while altitude(k) > 0
        % compute state
        a = -(2/tgo(:,k))*(Vf-V(:,k))...
            +(6/tgo(:,k)^2)*(rf-(r(:,k)+V(:,k)*tgo(:,k)));
        
        g = -mu*r(:,k)/((r(:,k)'*r(:,k))^(1.5));
        aT = a-g;
        
        % thrust limiter
        T_unlimited = norm(aT)*(m0-m_tot(k));
        
        if T_unlimited > T_max
            T(:,k) = T_max;
        elseif T_unlimited < T_min
            T(:,k) = T_min;
        else
            T(:,k) = T_unlimited;
        end
        limiter = T(:,k) / T_unlimited;
        
        aT = aT*limiter;
        
        
        % compute mass
        m_dot = T(:,k) / v_ex;
        m_tot(k+1) = m_tot(k) + m_dot*dt;
        
        % update tgo
        gamma = pi/2 - acos(r(:,k)'*V(:,k)/(norm(r(:,k))*speed(k)));
        gm = norm(g);
        
        a_tgo = 1/gm^2;
        b_tgo = sin(gamma)*speed(k)^2/(2*altitude(k)*gm^2);
        c_tgo = -(speed(k)^2*(1+sin(gamma)^2)/(4*altitude(k)*gm) + 1);
        
        a_GT = (-b_tgo + sqrt(b_tgo^2 - 4*a_tgo*c_tgo))/(2*a_tgo);
        if a_GT < 0
            a_GT = (-b_tgo - sqrt(b_tgo^2 - 4*a_tgo*c_tgo))/(2*a_tgo);
        end
        
        thrust_estimate(k) = norm(a_GT)*(m0-m_tot(k))*0.6/T_max;
        if thrust_estimate(k) >= 1
            thrust_estimate(k) = 1;
        elseif thrust_estimate(k) <= T_min/T_max
            thrust_estimate(k) = T_min/T_max;
        end
        
        p = 1000;
        tgo(k+1) = tgo(k);
        while p > tgo_tol
            tgo_temp = tgo(k+1);
            dVfV = Vf-V(:,k);
            dV = norm(dVfV + g*tgo_temp);
            tgo(k+1) = (m0/(thrust_estimate(k)*T_max/v_ex))*(1-exp(-dV/v_ex));
            p = norm(tgo(k+1) - tgo_temp);
        end
        % tgo(k+1) = 62-dt*k;
        
        % update state
        a = g + aT;
        V(:,k+1) = V(:,k) + a*dt;
        r(:,k+1) = r(:,k) + V(:,k)*dt;
        speed(k+1) = norm(V(:,k+1));
        altitude(k+1) = norm(r(:,k+1)) - R_M;
        range(k+1) = norm([r(2,k+1);r(3,k+1)]);
        
        k = k + 1;
    end
    
    % snip arrays to proper length
    r = r(:,1:k-1);
    V = V(:,1:k-1);
    T = T(1:k-1);
    tgo = tgo(1:k-1);
    m_tot = m_tot(1:k-1);
    speed = speed(1:k-1);
    altitude = altitude(1:k-1);
    range = range(1:k-1);
    thrust_estimate = thrust_estimate(1:k-1);
    
    %% Visualization
    
    trajectory = r-rf;
    
    % figure('Name','3D Trajectory')
    % plot3(trajectory(2,:),trajectory(3,:),trajectory(1,:))
    % xlabel('y')
    % ylabel('z')
    % zlabel('x')
    % title('Trajectory')
    %
    % figure('Name','tgo estimate')
    % plot(dt*[1:100:k-1],tgo(1:100:end))
    % title('t_go estimate')
    % xlabel('Run time (s)')
    %
    % figure('Name','Thrust')
    % plot(dt*[1:100:k-1],T(1:100:end)/T_max)
    % title('Commanded Thrust')
    % xlabel('Run time (s)')
    % ylabel('Throttle percentage')
    %
    % figure('Name','Fuel Consumption')
    % plot(dt*[1:100:k-1],m_tot(1:100:end))
    % title('Fuel Consumption')
    % xlabel('Run time (s)')
    % ylabel('Fuel Consumed')
    %
    % figure('Name','Speed')
    % plot(dt*[1:100:k-1],speed(1:100:end))
    % title('Speed')
    % xlabel('Run time (s)')
    %
    % figure('Name','Altitude')
    % plot(dt*[1:100:k-1],altitude(1:100:end))
    % title('Altitude')
    % xlabel('Run time (s)')
    %
    % figure('Name','Speed vs. Altitude')
    % plot(altitude(1:100:end),speed(1:100:end))
    % title('Speed vs. Altitude')
    % xlabel('Altitude')
    % ylabel('Speed')
    % set(gca,'Xdir','reverse')
    %
    % figure('Name','Range vs. Altitude')
    % plot(altitude(1:100:end),range(1:100:end))
    % title('Range vs. Altitude')
    % xlabel('Altitude')
    % ylabel('Range')
    % set(gca,'Xdir','reverse')
    
    fprintf('Case: %.0f\n', scenario)
    fprintf('Flight Time: %.2f (s)\n', (k-1)*dt)
    fprintf('Fuel: %.2f (kg)\n', m_tot(end))
    fprintf('Final range: %.2f (m)\n', range(end))
    fprintf('Final speed: %.2f (m/s)\n', speed(end))
    fprintf('Final angle to vertical: %.2f%c\n',...
        acosd(-aT'*g/(norm(aT)*norm(g))), char(176))
end

