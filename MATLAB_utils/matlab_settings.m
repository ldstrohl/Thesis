%% Matlab Settings
function matlab_settings()
% set print precision
format long g;
% default figures to docked
set(0,'DefaultFigureWindowStyle','docked');
% set color scheme
[~,result] = system('tasklist /FI "imagename eq matlab.exe"');
nMatlab = length(strfind(result,'.exe'));
switch nMatlab
    case 1
        setupSolarized('dark');
    case 2
        setupSolarized('light');
    otherwise
end
clear; 
end