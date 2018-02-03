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

% if(k>82250)
%     fprintf('k: %.0d\n',k)
%     fprintf('r: [%.2f %.2f %.2f]\n',r(1),r(2),r(3))
%     fprintf('V: [%.2f %.2f %.2f]\n',V(1),V(2),V(3))
%     fprintf('r_nav: [%.2f %.2f %.2f]\n',r_nav(1),r_nav(2),r_nav(3))
%     fprintf('V_nav: [%.2f %.2f %.2f]\n',V_nav(1),V_nav(2),V_nav(3))
%     fprintf('r_est: [%.2f %.2f %.2f]\n',r_est(1),r_est(2),r_est(3))
%     fprintf('V_est: [%.2f %.2f %.2f]\n\n',V_est(1),V_est(2),V_est(3))
% end



end

% Low Pass Filter
function x_est = LPF(x,x_prev)
alpha = 0.3;

x_est = alpha*x_prev + (1-alpha)*x;

end
