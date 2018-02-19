% Trajectory plot script
% Opens pre-generated 3-D Trajectory plot and overlays initial conditions
clear
close all
filename = 'traj_simvsadv.fig';
openfig(filename);

R_M = 3.39619*10^6;
hold on
for j = 6
    % Initial Condition
    switch j
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
%     af = af_factor*mu*rf/((rf'*rf)^(1.5));
    
    % coordinate transformation to landing site
    % [up;east;north]
    A = [cosd(lonf),sind(lonf),0;-sind(lonf),cosd(lonf),0;0,0,1];
    r0 = A*r0 - R_M*[1;0;0];
    V0 = A*V0;
    rf = A*rf - R_M*[1;0;0];
    Vf = A*Vf;
%     af = A*af;
    
r0 = [r0(2),r0(3),r0(1)];
rf = [rf(2),rf(3),rf(1)];
V0 = [V0(2),V0(3),V0(1)];
Vf = [Vf(2),Vf(3),Vf(1)];

plot3(r0(1),r0(2),r0(3),'X')

end
plot3(rf(1),rf(2),rf(3),'O')
grid on
h = gcf;

hold off