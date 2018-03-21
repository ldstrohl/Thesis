% Plot bullseye of final landing
clear
close all

%% inputs
savefig = 1;
figbase = 'bulldisppowatmo';

% data filename
rundata_file = 'rundata_atmo_bull.mat';

% define bullseye
inner_r = 50;
outer_r = 500;

%% initialization
load(rundata_file)

E = rundata(:,6);
N = rundata(:,7);
U = rundata(:,8);

% plot data
plot(rundata(:,6),rundata(:,7),'x')
hold on
circle([0 0],inner_r,'--')
circle([0 0],outer_r,'--')
xlabel('East (m)')
ylabel('North (m)')
legend('Final Landing Site',strcat(num2str(inner_r),' meter'),...
    strcat(num2str(outer_r),' meter'))
axis equal
grid on

if savefig
    thesis_fig(gcf,figbase)
end


function circle(center,r,style)
th = 0:pi/100:2*pi;
xunit = r * cos(th) + center(1);
yunit = r * sin(th) + center(2);
plot(xunit,yunit,style);
end
