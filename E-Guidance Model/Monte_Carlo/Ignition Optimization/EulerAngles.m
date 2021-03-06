function [pitch, yaw, roll, alpha, body_P] = EulerAngles(aT,r,V)
% Determines orientation of the spacecraft
% [up;east;north]

body_Y = aT/norm(aT);
body_P = cross(r,body_Y);
body_P = body_P/norm(body_P);
body_R = cross(body_Y,body_P);
pitch = asind(body_R(1));
yaw = atan2d(body_R(2),body_R(3));
roll = atan2d(body_P(1),body_Y(1));

alpha = atan2d(-(V'*body_Y),(V'*body_R));
end
