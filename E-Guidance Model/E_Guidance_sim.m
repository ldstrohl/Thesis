% E-Guidance simulation
% Lloyd Strohl
% 5/21/17

% Parameters
clear
close all
tgoest = 70;
R_M = 3.39619*10^6;

V0 = [82.5755978749391;...
    127.252297560431;...
    589.632069823669];
lonf = 184.2;   % azimuth
latf = 0;       % inclination
rf = [R_M*cosd(lonf)*cosd(latf);...
    R_M*sind(lonf)*cosd(latf);...
    R_M*sind(latf)]; 
Vf = (rf/norm(rf))*-1;
r0 = [-3394065.53680778;-253202.004401189;-19862.3977888315];
m0 = 58000;
T_max = 800000;
T_min = 200000;
v_ex = 3531.7;
mu = 4.282828185603917*10^13;
tol_tgo = 0.01;
rtol = 10;

% sim('a_T_vector_tgo')
sim('a_T_vector_tgo_iterator_v1')