clear
close all
alpha = [-90:1:90];
Mach = [0.2:.1:4];

cd = ones(length(alpha),length(Mach));
for i = 1:length(alpha)
    for j = 1:length(Mach)
        cd(i,j) = CD_supersonic(Mach(j),alpha(i));
    end
end

[X,Y] = meshgrid(Mach,alpha);
surf(X,Y,cd);
title('CD')
xlabel('Mach (deg)')
ylabel('\alpha (deg)')
zlabel('C_d')
        