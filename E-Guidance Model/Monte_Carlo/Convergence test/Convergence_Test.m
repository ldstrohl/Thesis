clear
close all

load('rundata_conv_RK4_thr.mat');
base = rundata(end,7);

for k = 1:length(rundata(:,1))-1
    dt(k) = rundata(k,end);
    err = (rundata(k,7) - base);
    RMSE(k) = norm(err);
end

y = dt.^1;
z = dt.^4*1e10;
loglog(dt,RMSE,dt,y, dt, z)
xlabel('\Delta t (s)')
ylabel('Speed Error (m/s)')
legend('\epsilon','O(\Delta t)', 'O(\Delta t^4)')
grid on
% thesis_fig(gcf,'convtestthr')