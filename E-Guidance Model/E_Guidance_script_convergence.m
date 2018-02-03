% E-Guidance - convergence test
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

% Rocket
m0 = 58000;
T_max = 800000;
T_min = 200000;
v_ex = 3531.7;


%% Simulation Settings
j = 1;
dtvec = [10^-1,10^-2,10^-3,10^-4];
for dt = dtvec % step size
    
    r_tol = 10^-1; % stopping criterion - distance from end point
    tgo_tol = 0.01;
    tgo0 = 100;
    
    % initialization
    r = r0*ones(1,tgo0/dt);
    V = V0*ones(1,tgo0/dt);
    T = zeros(1,tgo0/dt);
    a = zeros(3,1);
    aT = a;
    m_tot = zeros(1,tgo0/dt);
    tgo = tgo0*ones(1,tgo0/dt);
    
    t = 0;
    k = 2;
    
    %% Simulation loop
    while norm(rf-r(:,k-1)) > r_tol
        % compute state
        a = -(2/tgo(:,k))*(Vf-V(:,k))...
            +(6/tgo(:,k)^2)*(rf-(r(:,k)+V(:,k)*tgo(:,k)));
        V(:,k+1) = V(:,k) + a*dt;
        r(:,k+1) = r(:,k) + V(:,k)*dt;
        
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
        p = 1000;
        tgo(k+1) = tgo(k);
        while p > tgo_tol
            tgo_temp = tgo(k+1);
            dVfV = Vf-V(:,k);
            dV = norm(dVfV + g*tgo_temp);
            tgo(k+1) = (m0/m_dot)*(1-exp(-dV/v_ex));
            p = norm(tgo(k+1) - tgo_temp);
        end
        
        k = k + 1;
    end
    
    mtotiter(j) = m_tot(k-1);
    j = j + 1;
end

% % snip arrays to proper length
% r = r(:,1:k-1);
% V = V(:,1:k-1);
% T = T(1:k-1);
% tgo = tgo(1:k-1);
% m_tot = m_tot(1:k-1);
% 
% %% Visualization
% plot(dt*[1:10:k-1],V(1,1:10:end),...
%     dt*[1:10:k-1],V(2,1:10:end),...
%     dt*[1:10:k-1],V(3,1:10:end))
% 

loglog(dtvec,abs(mtotiter-mtotiter(end)))
xlabel('dt (s)');
ylabel('Mass Consumption Error')
title('Convergence')

poly = polyfit(log(dtvec(1:end-1)),...
    log(abs(mtotiter(1:end-1)-mtotiter(end))),1);
Convergence = poly(1)

