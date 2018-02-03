% E-Guidance
% fixed throttle tgo, optimization
% Lloyd Strohl
% 5/21/17

clear
close all

%% Conditions
% Mars
R_M = 3.39619*10^6;
mu = 4.282828185603917*10^13;

% Position
r0 = [-3394065.53680778;-253202.004401189;-19862.3977888315];
V0 = [82.5755978749391;...
    127.252297560431;...
    589.632069823669];

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

angle_iter = 1;
for theta = [linspace(-30,30,5)]
    Az = [cosd(theta),sind(theta),0;-sind(theta),cosd(theta),0;0,0,1];
    V0 = Az*V0;

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
fixed_thrust = [0.5, linspace(0.55,0.75,15), 0.9];

% initialization
r = r0*ones(1,tgo0/dt);
V = V0*ones(1,tgo0/dt);
T = ones(1,tgo0/dt);
speed = T;
altitude = T;
range = T;

a = zeros(3,1);
aT = a;
m_tot = zeros(1,tgo0/dt);
tgo = tgo0*ones(1,tgo0/dt);


    
    q = 1; % thrust fraction iterator
    for fixed_thrust_fraction = fixed_thrust
        %% Simulation loop
        t = 0;
        k = 2; % single simulation iterator
        while k < tgo0/dt && altitude(k-1) > 0
            % compute state
            a = -(2/tgo(:,k))*(Vf-V(:,k))...
                +(6/tgo(:,k)^2)*(rf-(r(:,k)+V(:,k)*tgo(:,k)));
            
            g = -mu*r(:,k)/((r(:,k)'*r(:,k))^(1.5));
            aT = a-g;
            speed(k) = norm(V(:,k));
            altitude(k) = norm(r(:,k)) - R_M;
            range(k) = norm([r(2,k);r(3,k)]);
            
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
            p = 1000;
            tgo(k+1) = tgo(k);
            while p > tgo_tol
                tgo_temp = tgo(k+1);
                dVfV = Vf-V(:,k);
                dV = norm(dVfV + g*tgo_temp);
                tgo(k+1) = (m0/(fixed_thrust_fraction*T_max/v_ex))*(1-exp(-dV/v_ex));
                p = norm(tgo(k+1) - tgo_temp);
            end
            %    tgo(k+1) = 62.1-dt*k;
            
            % update state
            a = g + aT;
            V(:,k+1) = V(:,k) + a*dt;
            r(:,k+1) = r(:,k) + V(:,k)*dt;
            
            k = k + 1;
        end
        
        % snip arrays to proper length
        r = r(:,2:k-1);
        V = V(:,2:k-1);
        T = T(2:k-1);
        tgo = tgo(2:k-1);
        m_tot = m_tot(2:k-1);
        speed = speed(2:k-1);
        altitude = altitude(2:k-1);
        range = range(2:k-1);
        
        % %% Visualization
        %
        % trajectory = r-rf;
        %
        % figure('Name','3D Trajectory')
        % plot3(trajectory(2,:),trajectory(3,:),trajectory(1,:))
        % xlabel('y')
        % ylabel('z')
        % zlabel('x')
        % title('Trajectory')
        %
        % figure('Name','tgo estimate')
        % plot(dt*[2:100:k-1],tgo(2:100:end))
        % title('t_go estimate')
        % xlabel('Run time (s)')
        %
        % figure('Name','Thrust')
        % plot(dt*[2:100:k-1],T(2:100:end)/T_max)
        % title('Commanded Thrust')
        % xlabel('Run time (s)')
        % ylabel('Throttle percentage')
        %
        % % figure('Name','Velocities')
        % % plot(dt*[1:10:k-1],V(1,1:10:end),...
        % %     dt*[1:10:k-1],V(2,1:10:end),...
        % %     dt*[1:10:k-1],V(3,1:10:end))
        % % title('Velocities')
        % % legend('x-velocity','y-velocity','z-velocity')
        %
        % % figure('Name','Coordinates')
        % % plot(dt*[1:10:k-1],r(1,1:10:end)-rf(1),...
        % %     dt*[1:10:k-1],r(2,1:10:end)-rf(2),...
        % %     dt*[1:10:k-1],r(3,1:10:end)-rf(3))
        % % title('Coordinates')
        % % legend('x','y','z')
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
        
        fraction(q) = fixed_thrust_fraction;
        flight_time(q) = (k-1)*dt;
        total_fuel(q) = m_tot(end);
        final_speed(q) = speed(end);
        
        %     fprintf('Thrust Fraction: %f\n', fraction(q))
        %     fprintf('Flight Time: %f\n', flight_time(q))
        %     fprintf('Fuel: %f\n', total_fuel(q))
        %     fprintf('Final speed: %f\n', final_speed(q))
        q = q + 1;
    end
    
    figure
    plot(fixed_thrust,total_fuel)
    title('Fuel Consumption vs. Average Throttle Estimate')
    xlabel('Throttle Estimate')
    ylabel('Fuel Consumed')
    
    total_fuel_variation(angle_iter) = fixed_thrust(total_fuel==min(total_fuel));
    angle_iter = angle_iter+1;
end
