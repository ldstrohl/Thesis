% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Low Pass Filter lpf(x)
% Lloyd Strohl
% 03/03/18
% Smooths discontinuities in signal x using a simple low pass filter
% inputs x = nx1 or 1xn vector 
% outputs y = nx1 or 1xn filtered vector
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = lpf(x)

% filter coefficient
c = 0.5;

y =  x;

for k = 2:length(x)
% new = old*(1-a)+cur*a
y(k) = x(k)*c + y(k-1)*(1-c);
end

    
    
end
