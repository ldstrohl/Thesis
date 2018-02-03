% E-Guidance
% Advanced - commanded final thrust acceleration
% Lloyd Strohl
% 5/21/17

clear
close all

%% Conditions
% Mars
R_M = 3.39619*10^6;
mu = 4.282828185603917*10^13;

for scenario = [1,2,3,4,5,6]
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
    
    lonf = 184.2;   % azimuth
    latf = 0;       % inclination
    rf = [R_M*cosd(lonf)*cosd(latf);...
        R_M*sind(lonf)*cosd(latf);...
        R_M*sind(latf)];
    Vf = (rf/norm(rf))*-1;
    af = mu*rf/((rf'*rf)^(1.5));
    
    % coordinate transformation to landing site
    A = [cosd(lonf),sind(lonf),0;-sind(lonf),cosd(lonf),0;0,0,1];
    r0 = A*r0;
    V0 = A*V0;
    rf = A*rf;
    Vf = A*Vf;
    af = A*af;
    
    % Euler Angles
    body_x = -V0/norm(V0);
    body_y = cross(body_x,r0);
    body_y = body_y/norm(body_y);
    body_z = cross(body_x,body_y);
    pitch = asind(body_x(1));
    yaw = atand(body_x(2)/body_x(3));
    roll = atand(body_y(1)/body_z(1));
    
    % Rocket parameters
    m0 = 58000;
    T_max = 800000;
    T_min = 200000;
    v_ex = 3531.7;
    
    % gravity turn tgo estimate
    gamma = pi/2 - acos(r0'*V0/(norm(r0)*norm(V0)));
    gm = norm(af);
    
    a_tgo = 1/gm^2;
    b_tgo = sin(gamma)*norm(V0)^2/(2*(norm(r0)-R_M)*gm^2);
    c_tgo = -(norm(V0)^2*(1+sin(gamma)^2)/(4*(norm(r0)-R_M)*gm) + 1);
    
    a_GT = (-b_tgo + sqrt(b_tgo^2 - 4*a_tgo*c_tgo))/(2*a_tgo);
    if a_GT < 0
        a_GT = (-b_tgo - sqrt(b_tgo^2 - 4*a_tgo*c_tgo))/(2*a_tgo);
    end
    
    tgo0 = round(norm(V0)/2 * ((1+sin(gamma))/(a_GT + gm)...
        + (1-sin(gamma))/(a_GT-gm)),0);
    
    % thrust_estimate = norm(a_GT)*(m0)*0.326/T_max;
    % % 1 - 3.0625/12581.98, 3.125/12513.88, 0.325/12460.10, 0.326/12459.62
    % % 2 - 3.0625/11604.10, 3.125/11449.34, 0.325/11214.66, 0.326/11199.05
    % % 3 - 3.0625/11651.23, 3.125/11444.26, 0.325/11167.98, 0.326/11151.05
    % % 4 - 3.0625/12465.41, 3.125/12068.73, 0.325/11627.03, 0.326/11604.68
    % % 5 - 3.0625/14011.14, 3.125/13197.29, 0.325/12376.66, 0.326/12341.77
    % % 6 - 3.0625/16199.89, 3.125/14836.89, 0.325/13320.87, 0.326/13272.87
    %
    % if thrust_estimate >= 1
    %     thrust_estimate = 1;
    % elseif thrust_estimate <= T_min/T_max
    %     thrust_estimate = T_min/T_max;
    % end
    
    %% Simulation Settings
    dt = 10^-3; % step size
    tgo_tol = 0.01;
    tgo_stop = 0.5;
    throttle_estimate = .305;
    
    % initialization
    r = r0*ones(1,tgo0/dt);
    V = V0*ones(1,tgo0/dt);
    T = ones(1,tgo0/dt);
    speed = T*norm(V(:,1));
    altitude = T*norm(r(:,1)) - R_M;
    range = T*norm([r(2,1);r(3,1)]);
    pitch = pitch*T;
    yaw = yaw*T;
    roll = roll*T;
    alpha = T;
    
    a = zeros(3,1);
    aT = a;
    m_tot = zeros(1,tgo0/dt);
    tgo = tgo0*ones(1,tgo0/dt);
    
    t = 0;
    k = 1;
    
    %% Simulation loop
    while altitude(k) > 0
        % compute acceleration
        g = -mu*r(:,k)/((r(:,k)'*r(:,k))^(1.5));
        
        Ak = [0,0,1;...
            18/tgo(k)^2,-24/tgo(k)^3,-6/tgo(k);...
            -24/tgo(k)^3,36/tgo(k)^4,6/tgo(k)^2];
        
        % x-constants
        FVx = [Vf(1)-V(1,k);...
            rf(1)-r(1,k)-V(1,k)*tgo(k);...
            g(1)+af(1)];
        kx = Ak*FVx;
        % y-constants
        FVy = [Vf(2)-V(2,k);...
            rf(2)-r(2,k)-V(2,k)*tgo(k);...
            g(2)+af(2)];
        ky = Ak*FVy;
        
        % z-constants
        FVz = [Vf(3)-V(3,k);...
            rf(3)-r(3,k)-V(3,k)*tgo(k);...
            g(3)+af(3)];
        kz = Ak*FVz;
        
        k1 = [kx(1);ky(1);kz(1)];
        k2 = [kx(2);ky(2);kz(2)];
        k3 = [kx(3);ky(3);kz(3)];
        
        a = k1+k2*tgo(k) + k3*tgo(k)^2;
        
        
       if tgo(k) > tgo_stop % hold aT constant near end
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
        else 
            T_unlimited = norm(aT)*(m0-m_tot(k));
            if T_unlimited > T_max
                T(:,k) = T_max;
            elseif T_unlimited < T_min
                T(:,k) = T_min;
            else
                T(:,k) = T_unlimited;
            end
            limiter = T(:,k) / T_unlimited;
        end
        
        
        % compute mass
        m_dot = T(:,k) / v_ex;
        m_tot(k+1) = m_tot(k) + m_dot*dt;
        
        % update tgo
        p = 1000;
        tgo(k+1) = tgo(k);
        iter = 0;
        while p > tgo_tol && iter < 100
            tgo_temp = tgo(k+1);
            dVfV = Vf-V(:,k);
            dV = norm(dVfV + g*tgo_temp);
            tgo(k+1) = (m0/(throttle_estimate*T_max/v_ex))*(1-exp(-dV/v_ex));
            p = norm(tgo(k+1) - tgo_temp);
            iter = iter + 1;
        end
        %    tgo(k+1) = 62.1-dt*k;
        
        % update state
        body_x = aT/norm(aT);
        body_y = cross(body_x,r(:,k));
        body_y = body_y/norm(body_y);
        body_z = cross(body_x,body_y);
        pitch(k) = asind(body_x(1));
        yaw(k) = atand(body_x(2)/body_x(3));
        roll(k) = atand(body_y(1)/body_z(1));
        alpha(k) = atand((V(:,k)'*body_z)/(V(:,k)'*body_x));
        
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
    pitch = pitch(1:k-1);
    yaw = yaw(1:k-1);
    roll = roll(1:k-1);
    alpha = alpha(1:k-1);
    
    %% Visualization
    
    trajectory = r-rf;
    
    figure('Name','3D Trajectory')
    plot3(trajectory(2,:),trajectory(3,:),trajectory(1,:))
    xlabel('y')
    ylabel('z')
    zlabel('x')
    title('Trajectory')
    
    figure('Name','tgo estimate')
    plot(dt*[1:100:k-1],tgo(1:100:end))
    title('t_go estimate')
    xlabel('Run time (s)')
    
    figure('Name','Thrust')
    plot(dt*[1:100:k-1],T(1:100:end)/T_max)
    title('Commanded Thrust')
    xlabel('Run time (s)')
    ylabel('Throttle percentage')
    
    figure('Name','Fuel Consumption')
    plot(dt*[1:100:k-1],m_tot(1:100:end))
    title('Fuel Consumption')
    xlabel('Run time (s)')
    ylabel('Fuel Consumed')
    
    figure('Name','Speed')
    plot(dt*[1:k-1],speed)
    title('Speed')
    xlabel('Run time (s)')
    
    figure('Name','Altitude')
    plot(dt*[1:100:k-1],altitude(1:100:end))
    title('Altitude')
    xlabel('Run time (s)')
    
    figure('Name','Speed vs. Altitude')
    plot(altitude(1:100:end),speed(1:100:end))
    title('Speed vs. Altitude')
    xlabel('Altitude')
    ylabel('Speed')
    set(gca,'Xdir','reverse')
    
    figure('Name','Range vs. Altitude')
    plot(altitude(1:100:end),range(1:100:end))
    title('Range vs. Altitude')
    xlabel('Altitude')
    ylabel('Range')
    set(gca,'Xdir','reverse')
    figure('Name','Pitch')
    plot(dt*[1:k-1],pitch)
    title('Pitch')
    xlabel('Run time (s)')
    
    %% Performance messages
    fprintf('Case: %.0f\n', scenario)
    fprintf('Flight Time: %.2f (s)\n', (k-1)*dt)
    fprintf('Fuel: %.2f (kg)\n', m_tot(end))
    fprintf('Final range: %.2f (m)\n', range(end))
    fprintf('Final speed: %.2f (m/s)\n', speed(end))
    fprintf('Final angle to vertical: %.2f%c\n',...
        acosd(-aT'*g/(norm(aT)*norm(g))), char(176))
    % fprintf('Final Thrust x-cord %.2f\n',aT(1))
    
end