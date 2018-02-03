% Integration Module
% inputs: current state, dt
% outputs: next state

function [r1,V1] = integration_module(r0,V0,a_TLD,mu,dt)
% g = -mu*r0/((r0'*r0)^(1.5));
% a0 = g+aT;
% 
% % first order explicit Euler 
% V1 = V0 + a0*dt;
% r1 = r0 + V0*dt;

% RK4
phi0 = [r0;V0];

k1 = f(phi0,a_TLD,mu);
k2 = f(phi0 + (dt/2)*k1,a_TLD,mu);
k3 = f(phi0 + (dt/2)*k2,a_TLD,mu);
k4 = f(phi0 + dt*k3,a_TLD,mu);

% % error checking
% e1 = dt*sqrt(k1.*k1);
% if any(e1>2)
%     for j = 1:6
%         fprintf('e1(%.0f): %.2f\n',j,e1(j))
%     end
%     fprintf('|F|: %.2f\n\n',norm(F))
% end
% % e2 = dt*sqrt(k2.*k2);
% % e3 = dt*sqrt(k3.*k3);
% % e4 = dt*sqrt(k4.*k4);


    


phi1 = phi0 + (dt/6)*(k1+2*k2+2*k3+k4);

r1 = phi1(1:3);
V1 = phi1(4:6);

end

function dphidt = f(phi,a_TLD,mu)
r = phi(1:3);
V = phi(4:6);
g = -mu*r/((r'*r)^(1.5));

dphidt = [V;g+a_TLD];
end
