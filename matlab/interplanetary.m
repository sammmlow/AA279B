% Physical Constants for Problem 8.3
% AA 279A Winter 2022
% Andrew K. Barrows

clc; clear all; close all;

% Physical constants from https://ssd.jpl.nasa.gov/horizons/
% Valid for 0000h 1 March 2022 CT

% Sun initial conditions
musun = 132712440041.93938; % [km^3/sec^2]
rsun0 = [ 1.036867385628283E+06; ...
         -3.164409455542757E+05; ...
         -2.755319914038984E+04];
vsun0 = [ 2.577398212280307E-03; ...
          1.450959835015759E-02; ...
         -1.406783318786790E-04];

% Venus initial conditions
muvenus = 324858.592; % [km^3/sec^2]
rvenus0 = [ 3.486376123235895E+06; ...
           -6.925254981998456E+07; ...
           -5.887103744782358E+06];       
vvenus0 = [ 3.891394101959566E+01; ...
            4.239420589378792E+00; ...
           -3.222229903892365E+00];

% % Earth initial conditions
muearth = 398600.435436; % [km^3/sec^2] 
rearth0 = [-2.351483767168928E+07; ...
            1.447216737998444E+08; ...
           -4.073146757017076E+04];  
vearth0 = [-2.984201918638308E+01; ...
           -5.064455455220450E+00; ...
           -1.614822049471609E-04];

% Moon initial conditions
mumoon = 4902.800066; % [km^3/sec^2]
rmoon0 = [-2.368927882124358E+07; ...
           1.450859584686864E+08; ...
          -1.543499608439207E+04];
vmoon0 = [-3.071086386096760E+01; ...
          -5.503343770995965E+00; ...
           6.331719512134981E-02];

% Saturn initial conditions
musat = 37931206.234; % [km^3/sec^2]
rsat0 = [-1.419649252421490E+09; ...
         -1.292355436624920E+08; ...
          5.876654433717017E+07];
vsat0 = [ 3.501943867478344E-01; ...
         -9.636829713639742E+00; ...
          1.530884311611072E-01];

% Jupiter initial conditions
mujup = 126686531.900; % [km^3/sec^2]
rjup0 = [-7.999722358881212E+08; ...
          1.324911728292180E+08; ...
          1.734264413707683E+07];
vjup0 = [-2.294859734285312E+00; ...
         -1.228255093370500E+01; ...
          1.024347284634670E-01];



%% Main code below.
% Initialize state
% y0 = [ musun rsun0' vsun0' ...
%        mumerc rmerc0' vmerc0' ...
%        muvenus rvenus0' vvenus0' ...
%        muearth rearth0' vearth0' ...
%        mumoon rmoon0' vmoon0' ];


y0 = [ musun rsun0' vsun0' ...
       muvenus rvenus0' vvenus0' ...
       muearth rearth0' vearth0' ...
       mumoon rmoon0' vmoon0' ...
       musat rsat0' vsat0' ...
       mujup rjup0' vjup0' ];

% z = nbody_statedot( [0], y0 ) % Test line

options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
[t, y] = ode113(@nbody_statedot, [0:3600:315360000]', y0, options);

% Plot the line trajectory
figure(1);
plot( y(:, 2), y(:, 3), LineWidth=2 );                 % Sun
hold('on'); grid('on');
plot( y(:, 9), y(:,10), LineWidth=2 );                 % Venus
plot( y(:,16), y(:,17), LineWidth=2 );                 % Earth
plot( y(:,23), y(:,24), LineWidth=2, LineStyle='--' ); % Moon
plot( y(:,30), y(:,31), LineWidth=2, LineStyle='--' ); % Saturn
plot( y(:,37), y(:,38), LineWidth=2, LineStyle='--' ); % Jupiter
xlabel('X-Coordinate [km]')
ylabel('Y-Coordinate [km]')
title('Plot of Inner Solar System Bodies for 1 Year')

% Plot the end-points of their trajectory
scatter( y(end, 2), y(end, 3), 'k', LineWidth=1.5 ); % Sun
scatter( y(end, 9), y(end,10), 'k', LineWidth=1.5 ); % Venus
scatter( y(end,16), y(end,17), 'k', LineWidth=1.5 ); % Earth
scatter( y(end,23), y(end,24), 'k', LineWidth=1.5 ); % Moon
scatter( y(end,30), y(end,31), 'k', LineWidth=1.5 ); % Saturn
scatter( y(end,37), y(end,38), 'k', LineWidth=1.5 ); % Jupiter


legend('Sun','Venus','Earth','Moon', ...
       'Saturn','Jupiter','','','','','','');

% % Animate trajectory
% hold('on'); 
% comet( y(:, 2), y(:, 3) ); % Sun
% hold('on'); 
% comet( y(:, 9), y(:,10) ); % Venus
% comet( y(:,16), y(:,17) ); % Earth
% comet( y(:,23), y(:,24) ); % Moon
% comet( y(:,30), y(:,31) ); % Saturn
% comet( y(:,37), y(:,38) ); % Jupiter



%% Function that returns the state dot of the entire solar system
% state = [ musun rsun0' vsun0' ...
%           mumerc rmerc0' vmerc0' ...
%           muvenus rvenus0' vvenus0' ...
%           muearth rearth0' vearth0' ...
%           mumoon rmoon0' vmoon0' ]';
function [statedot] = nbody_statedot( t, state )
    L = length(state);     % Number of states
    statedot = zeros(L,1); % Output state
    for k = 1:7:L
        statedot(k+1) = state(k+4);
        statedot(k+2) = state(k+5);
        statedot(k+3) = state(k+6);
        accel = zeros(3,1);
        for j = 1:7:L
            if j ~= k
                mu = state(j);
                r = state(k+1:k+3) - state(j+1:j+3);
                accel = accel - (mu * r / (norm(r)^3));
            end
        end
        statedot(k+4:k+6) = accel;
    end
end