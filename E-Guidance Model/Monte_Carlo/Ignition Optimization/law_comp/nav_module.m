% navigation module
% inputs: states
% output: states plus noise
% Module takes current state from simulation and generates zero mean
% Gaussian noise, adds them, and sends them back as navigation state
% for use in the guidance law

function [r_est,V_est] = nav_module(r,V,r_est_prev,V_est_prev,r_sig,V_sig)

r_nav = r + randn(1)*r_sig;
V_nav = V + randn(1)*V_sig;

r_est = LPF(r_nav,r_est_prev);
V_est = LPF(V_nav,V_est_prev);

end

% Low Pass Filter
function x_est = LPF(x,x_prev)
alpha = 0.3;

x_est = alpha*x_prev + (1-alpha)*x;

end
