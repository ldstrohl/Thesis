% E-Guidance
% Souza tgo
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


% Rocket parameters
m0 = 58000;
T_max = 800000;
T_min = 200000;
v_ex = 3531.7;


%% Simulation Settings
dt = 10^-3; % step size
tgo0 = 200;

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

t = 0;
k = 2; 

%% Simulation loop
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
   
   % update state
   a = g + aT;
   V(:,k+1) = V(:,k) + a*dt;
   r(:,k+1) = r(:,k) + V(:,k)*dt;
   
%    % update tgo
   p = [1.5*10^7+norm(g)^2/2,...
       0,...
       -2*(V(:,k)'*V(:,k)),...
       -1,...
       -12*V(:,k)'*r(:,k) - 18*(r(:,k)'*r(:,k))];
   tgo_roots = roots(p);
   x = ones(1,length(tgo_roots));
   for j = 1:length(tgo_roots)
       x(j) = isreal(tgo_roots(j));
   end
   p = tgo_roots(x>0);
   tgo(k+1) = p(p>0);


%    a_tgo = -2*(V(:,k)'*V(:,k))/(norm(g)^2/2);
%    b_tgo = -12*(V(:,k)'*(-r(:,k)+rf))/(norm(g)^2/2);
%    c_tgo = -18*((-r(:,k)+rf)'*(-r(:,k)+rf))/(norm(g)^2/2);
%    beta = (1/27)*(16*a_tgo^3 - 18*a_tgo*(a_tgo^2-4*c_tgo)-27*b_tgo^2);
%    alpha = (1/3)*(3*(a_tgo^2-4*c_tgo)-4*a_tgo^2);
%    delta = alpha^3/27 + beta^2/4;
%    Z = nthroot(-(beta/2) + sqrt(delta),3) + nthroot(-(beta/2) - sqrt(delta),3);
%    eta = Z - 2*a_tgo/3;
%    zeta = -b_tgo/(2*eta);
%    xi = (a_tgo + eta)/2;
%    
%    tgo_sqrt_value = eta - 4*(xi-sqrt(eta)*zeta);
%    if tgo_sqrt_value < 0
%        tgo_sqrt_value = eta - 4*(xi+sqrt(eta)*zeta);
%        tgo_numerator = -sqrt(eta) + sqrt(tgo_sqrt_value);
%        if tgo_numerator < 0
%            tgo(k+1) = (-sqrt(eta) - sqrt(tgo_sqrt_value))/2;
%        else
%            tgo(k+1) = tgo_numerator/2;
%        end
%    else
%        tgo(k+1) = (sqrt(eta)+sqrt(tgo_sqrt_value))/2;
%        if tgo(k+1) < 0
%            tgo(k+1) = (sqrt(eta)+sqrt(tgo_sqrt_value))/2;
%        end
%    end

% tgo1 = (sqrt(eta) + sqrt(eta - 4*(xi-sqrt(eta)*zeta)))/2;
% tgo2 = (sqrt(eta) - sqrt(eta - 4*(xi-sqrt(eta)*zeta)))/2;   
   
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

%% Visualization

trajectory = r-rf;

figure('Name','3D Trajectory')
plot3(trajectory(2,:),trajectory(3,:),trajectory(1,:))
xlabel('y')
ylabel('z')
zlabel('x')
title('Trajectory')

figure('Name','tgo estimate')
plot(dt*[2:100:k-1],tgo(2:100:end))
title('t_go estimate')
xlabel('Run time (s)')

figure('Name','Thrust')
plot(dt*[2:100:k-1],T(2:100:end))
title('Commanded Thrust')
xlabel('Run time (s)')

figure('Name','Velocities')
plot(dt*[1:10:k-1],V(1,1:10:end),...
    dt*[1:10:k-1],V(2,1:10:end),...
    dt*[1:10:k-1],V(3,1:10:end))
title('Velocities')
legend('x-velocity','y-velocity','z-velocity')

% figure('Name','Coordinates')
% plot(dt*[1:10:k-1],r(1,1:10:end)-rf(1),...
%     dt*[1:10:k-1],r(2,1:10:end)-rf(2),...
%     dt*[1:10:k-1],r(3,1:10:end)-rf(3))
% title('Coordinates')
% legend('x','y','z')

figure('Name','Fuel Consumption')
plot(dt*[1:100:k-1],m_tot(1:100:end))
title('Fuel Consumption')
xlabel('Run time (s)')
ylabel('Fuel Consumed')

figure('Name','Speed')
plot(dt*[1:100:k-1],speed(1:100:end))
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

fprintf('Flight Time: %.2f\n', (k-1)*dt)
fprintf('Fuel: %.2f\n', m_tot(end))
fprintf('Final range: %.2f\n', range(end))
fprintf('Final speed: %.2f\n', speed(end))