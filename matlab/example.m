clc; clear all; close all;
set(groot, 'defaultLineLineWidth', 2);
addpath("library");

%% Working example for a single call for HCI coords of Saturn and Titan
AU = 149597870.7; % [km]
epoch = datetime(2040, 6, 15, 12, 0, 0);
pv_saturn = ephemeris( 'S', epoch );
pv_titan  = ephemeris( 'T', epoch );

%% Example plot of trajectories of planets and moons for 5 years.
% 06x Objects tracked: Earth, Moon, Mars, Jupiter, Saturn, Titan
bodies      = 6;
timestep    = 7200;
timescale   = 0 : timestep : (5*365*86400);
timesamples = length(timescale);
eph         = zeros(6, timesamples, bodies);

for k = 1 : timesamples
    eph(:,k,1) = ephemeris( 'E', epoch ); % Earth
    eph(:,k,2) = ephemeris( 'L', epoch ); % Moon
    eph(:,k,3) = ephemeris( 'M', epoch ); % Mars
    eph(:,k,4) = ephemeris( 'J', epoch ); % Jupiter
    eph(:,k,5) = ephemeris( 'S', epoch ); % Saturn
    eph(:,k,6) = ephemeris( 'T', epoch ); % Titan
    epoch = epoch + seconds(timestep);
end

%% Plot the motion of all celestial bodies!
figure();
scatter([0], [0], 48, 'yellow', 'filled');
hold on; grid on;
for b = 1 : bodies
    plot( eph(2,:,b)./AU, eph(1,:,b)./AU );
    scatter( eph(2,end,b)./AU, eph(1,end,b)./AU, 'k', 'filled');
end
ylabel('X [Heliocentric, AU]');
xlabel('Y [Heliocentric, AU]');
set(gca,'XDir','reverse');
axis equal;


%% Test a random Lambert-based trajectory solver

trajectory = {'E','J','S'};
test_dates = [datetime(2030,1,1,0,0,0), ...
              datetime(2035,1,1,0,0,0), ...
              datetime(2040,1,1,0,0,0)];

results = computeTrajectory( trajectory, test_dates );