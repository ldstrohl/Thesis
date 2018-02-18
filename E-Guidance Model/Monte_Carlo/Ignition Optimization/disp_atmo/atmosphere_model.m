function [L,D, C_drag] = atmosphere_model(alt,V,alpha,S,body_P)
% computes Lift and Drag vectors in Mars atmosphere
% Inputs:   alt = altitude (m)
%           V = velocity (m/s)
%           alpha = angle of attack (deg)
%           body_P = body pitch axis vector
normV = norm(V);

% get Mars properties
[rho,~] = density(alt/1000);
Vs = Vsound(alt/1000);
Mach = normV/Vs;

% Compute lift and drag
C_lift = CL(Mach,alpha);
% C_lift = 1;
% C_drag = CD(Mach,alpha);
C_drag = CD(Mach,alpha);

% C_drag = 2;
L = cross(body_P,V);
L = L / norm(L);
D = -V / normV;
L = L * (1/2) * rho * normV^2 * S * C_lift;
D = D * (1/2) * rho * normV^2 * S * C_drag;


end
