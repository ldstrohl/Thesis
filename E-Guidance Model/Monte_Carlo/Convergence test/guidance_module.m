% Guidance Module
% Inputs: state, final condition, tgo, law option
% outputs: commanded throttle in direction of commanded thrust vector

function throttle = guidance_module(r,V,g,rf,Vf,tgo,m,T_max,law,af)

switch law
    case 1
        aT = -(2/tgo)*(Vf-V)...
            +(6/tgo^2)*(rf-(r+V*tgo)) - g;
    case 2
        Ak = [0,0,1;...
            18/tgo^2,-24/tgo^3,-6/tgo;...
            -24/tgo^3,36/tgo^4,6/tgo^2];
        
        % x-constants
        FVx = [Vf(1)-V(1);...
            rf(1)-r(1)-V(1)*tgo;...
            g(1)+af(1)];
        kx = Ak*FVx;
        % y-constants
        FVy = [Vf(2)-V(2);...
            rf(2)-r(2)-V(2)*tgo;...
            g(2)+af(2)];
        ky = Ak*FVy;
        
        % z-constants
        FVz = [Vf(3)-V(3);...
            rf(3)-r(3)-V(3)*tgo;...
            g(3)+af(3)];
        kz = Ak*FVz;
        
        k1 = [kx(1);ky(1);kz(1)];
        k2 = [kx(2);ky(2);kz(2)];
        k3 = [kx(3);ky(3);kz(3)];
        
        aT = k1+k2*tgo + k3*tgo^2 - g;
end
throttle = aT*m/T_max;

end
