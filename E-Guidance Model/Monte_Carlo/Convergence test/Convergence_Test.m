clear
close all

load('rundata_conv_RK4.mat');
base = rundata(end,7);

for k = 17:length(rundata(:,1))-1
    dt(k) = rundata(k,end);
    err = (rundata(k,7) - base);
    RMSE(k) = norm(err);
end

y = dt.^2*1e6;
loglog(dt,RMSE,dt,y)
xlabel('\Delta t (s)')
ylabel('Fuel Error (kg)')
legend('\epsilon','O(\Delta t)')