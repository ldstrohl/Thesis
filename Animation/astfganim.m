%% Create a Flight Animation from Trajectory Data
% This example shows how to create a flight animation for a trajectory using
% a FlightGear Animation object.
%
% *Note:* When running this example within the product, you must customize 
% the example with your FlightGear installation and uncomment the
% GenerateRunScript, system and play commands. You must also copy the
% $MATLAB/toolbox/aero/astdemos/HL20 folder into the
% $FLIGHTGEAR/data/Aircraft/ folder.  

% Copyright 1990-2017 The MathWorks, Inc.

%% Load Recorded Flight Trajectory Data
%
% The flight trajectory data for this example is stored in a comma
% separated value formatted file.  Use *dlmread* to read the data from 
% the file starting at row 1 and column 0 skipping the header information.

% [time lat lon alt 
% tdata = dlmread('asthl20log.csv',',',1,0);
% tdata(:,4) = tdata(:,4)+ones(length(tdata(:,4)),1)*20;
%% Create a Time Series Object from Trajectory Data
%
% Use the MATLAB(R) *timeseries* command to create the time series object, 
% ts, from the latitude, longitude, altitude, and Euler angle data along 
% with the time array in tdata.  To convert the latitude, longitude, and 
% Euler angles from degrees to radians use the *convang* function. 

ts = timeseries([convang(tdata(:,[3 2]),'deg','rad') ...
                 tdata(:,4) convang(tdata(:,5:7),'deg','rad')],tdata(:,1));

%%
% You can create imported data from this data using other valid formats,
% such as 'Array6DoF', For example:
%
% ts = [tdata(:,1) convang(tdata(:,[3 2]),'deg','rad') tdata(:,4) ...
%                                        convang(tdata(:,5:7),'deg','rad')];
%
% and 'Array3DoF'.
%
% ts = [tdata(:,1) convang(tdata(:,3),'deg','rad') tdata(:,4) ...
%                                         convang(tdata(:,6),'deg','rad')];

%% Use FlightGearAnimation Object to Initialize Flight Animation

%%
% Open a FlightGearAnimation object.
%%
h = Aero.FlightGearAnimation;

%%
% Set FlightGearAnimation object properties for timeseries.
%%

h.TimeseriesSourceType = 'Timeseries';
h.TimeseriesSource = ts;

%%
% Set FlightGearAnimation object properties about FlightGear
%%
% These properties include the path to the installation folder, the version 
% number, the aircraft geometry model, and the network information for
% FlightGear flight simulator. 

h.FlightGearBaseDirectory = 'C:\Program Files\FlightGear';
h.FlightGearVersion = '2016.1';
h.GeometryModelName = 'HL20';
h.DestinationIpAddress = '127.0.0.1';
h.DestinationPort = '5502';

%%
% Set the desired initial conditions (location and orientation) for
% FlightGear flight simulator.  

h.AirportId = 'KSFO';
h.RunwayId = '10L';
h.InitialAltitude = 7224;
h.InitialHeading = 113;
h.OffsetDistance = 4.72;
h.OffsetAzimuth = 0;

%%
% Set the seconds of animation data per second of wall-clock time.
h.TimeScaling = 10;

%%
% Use get(h) to check the FlightGearAnimation object properties and their
% values.
get(h)

%% Create a Run Script to Launch FlightGear Flight Simulator
%%
% To start FlightGear with the desired initial conditions (location, date,
% time, weather, and operating modes), create a run script with
% the *GenerateRunScript* command.  By default, *GenerateRunScript*
% saves the run script as a text file named 'runfg.bat'.

%%
% GenerateRunScript(h)

%%
% You do not need to generate this file each time the data is viewed.
% Generate it only when the desired initial conditions or FlightGear information
% changes.

%% Start FlightGear Flight Simulator
%% 
% To start FlightGear from the MATLAB command prompt, type the *system* command 
% to execute the run script created by *GenerateRunScript*.

%%
% system('runfg.bat &');

%%
% *Tip:* With the FlightGear window in focus, press the V key to alternate
% between the different aircraft views: cockpit view, helicopter view, and
% chase view. 

%% Play the Flight Animation of Trajectory Data
%%
% Once FlightGear is up and running, the FlightGearAnimation object can
% start to communicate with FlightGear.  To display the flight animation 
% with FlightGear, use the *play* command.

%%
% play(h)

%%
% % To display a screenshot of the flight animation, use the MATLAB *image* command.
% image(imread([matlabroot filesep fullfile('toolbox','aero','astdemos','figures','astfganim01.png')],'png'));
% axis off;
% set(gca,'Position',[ 0 0 1 1 ]);
% set(gcf,'MenuBar','none');
% 
% displayEndOfDemoMessage(mfilename)
