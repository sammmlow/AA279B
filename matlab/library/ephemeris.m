function hci = ephemeris( body, epoch )
% Inputs:
%   - body  --> char array e.g. 'Jupiter'
%   - epoch --> datetime object e.g. datetime(2000,1,1,12,0,0)
% Output:
%   - hci   --> 6x1 heliocentric inertial coords [km, km/s]

% Ephemeris source: http://ssd.jpl.nasa.gov/?planet_pos
% Validity: 1800 - 2050 AD

ID = 0;
AU = 149597870.7;          % [km]
GM_SUN = 1.32712440e11;    % [km^3/s^2]
GM_PLANET = nan;           % [km^3/s^2]

epoch_ref = datetime(2000, 1, 1, 12, 0, 0); % J2000 epoch
dt_days = days(epoch - epoch_ref);
dt_seconds = dt_days * 86400;

% ============================================================
% COMPUTATION OF PLANET EPHEMERIS IN HCI COORDINATES
% ============================================================

% Tag the body with an ID. Designed to be non-case sensitive.
switch body
    case {'Me'}      % Mercury
        ID = 1; GM_PLANET = 2.2032e4;
    case {'V'}       % Venus
        ID = 2; GM_PLANET = 3.24859e5;
    case {'E', 'L'}  % Earth or the Moon
        ID = 3; GM_PLANET = 3.98600e5;
    case {'M'}       % Mars
        ID = 4; GM_PLANET = 4.28283e4;
    case {'J'}       % Jupiter
        ID = 5; GM_PLANET = 1.26686534e10;
    case {'S', 'T'}  % Saturn or Titan
        ID = 6; GM_PLANET = 3.7931187e9;
    case {'U'}       % Uranus
        ID = 7; GM_PLANET = 5.793939e8;
    case {'N'}       % Neptune
        ID = 8; GM_PLANET = 6.836529e8;
    case {'P'}       % Pluto
        ID = 9; GM_PLANET = 8.719e2;
    otherwise
        warning('Input body is unexpected! Exiting');
        pos = nan(3,1);
        vel = nan(3,1);
        hci = [pos ; vel];
        return;
end

% Heliocentric orbital elements of planets for ID 1 to 9. 
% a = semimajor axis [AU]
% e = eccentricity
% i = inclination [deg]
% O = longitude of the ascending node [deg]
% o = longitude of periapsis [deg]
% L = mean longitude [deg]

planet_ephemeris = ...
[0.38709927 0.20563593  7.00497902  48.33076593  77.45779628 252.25032350; 
 0.72333566 0.00677672  3.39467605  76.67984255 131.60246718 181.97909950;  
 1.00000261 0.01671123 -0.00001531   0.00000000 102.93768193 100.46457166;  
 1.52371034 0.09339410  1.84969142  49.55953891 -23.94362959  -4.55343205; 
 5.20288700 0.04838624  1.30439695 100.47390909  14.72847983  34.39644501; 
 9.53667594 0.05386179  2.48599187 113.66242448  92.59887831  49.95424423; 
19.18916464 0.04725744  0.77263783  74.01692503 170.95427630 313.23810451; 
30.06992276 0.00859048  1.77004347 131.78422574  44.96476227 -55.12002969;  
39.48211675 0.24882730 17.14001206 110.30393684 224.06891629 238.92903833];

% Get specific heliocentric orbital elements of the planet.
planet_elements = planet_ephemeris(ID, :);

a = planet_elements(1) * AU;
e = planet_elements(2);
i = deg2rad(planet_elements(3));
O = deg2rad(planet_elements(4));
w = deg2rad(planet_elements(5) - planet_elements(4));
M = deg2rad(planet_elements(6) - planet_elements(5));

% Propagate the two-body motion of the planet.
n = sqrt(GM_SUN/(a^3));

% Propagated mean, eccentric, and true anomalies.
Mf = wrapTo2Pi(M + dt_seconds * n); 
Ef = M2E( Mf, e, 1e-6 );
vf = E2T( Ef, e );
planet_elements_f = [ a, e, i, O, w, vf ]';

% Compute heliocentric inertial (HCI) coordinates
[pos, vel] = OE2XCI( planet_elements_f, GM_SUN );
hci = [pos ; vel];

% ============================================================
% COMPUTATION OF PLANETARY MOONS' EPHEMERIS IN HCI COORDINATES
% ============================================================

ID_MOON = 0;

switch body
    case {'L'}  % Moon
        ID_MOON = 1;
    case {'T'}  % Titan
        ID_MOON = 2;
    otherwise
        return;
end

% Planetocentric orbital elements of moons. 
% aS = semimajor axis [AU]
% eS = eccentricity
% iS = inclination [deg]
% OS = longitude of the ascending node [deg]
% oS = argument of periapsis [deg]
% MS = mean anomaly [deg]

satellite_ephemeris = ...
[0.002569555 0.05540 5.16000 125.08 318.15 135.27 ;
 0.008167897 0.02920 0.34854  78.60  78.30  11.70 ];

% Get specific planetocentric orbital elements of the planet.
satellite_elements = satellite_ephemeris(ID_MOON, :);

aS = satellite_elements(1) * AU;
eS = satellite_elements(2);
iS = deg2rad(satellite_elements(3));
OS = deg2rad(satellite_elements(4));
wS = deg2rad(satellite_elements(5));
MS = deg2rad(satellite_elements(6));

% Propagate the two-body motion of the planet.
nS = sqrt(GM_PLANET/(aS^3));

% Propagated mean, eccentric, and true anomalies.
MSf = wrapTo2Pi(MS + dt_seconds * nS);
ESf = M2E( MSf, eS, 1e-6 );
vSf = E2T( ESf, eS );
satellite_elements_f = [ aS, eS, iS, OS, wS, vSf ]';

% Compute heliocentric inertial (HCI) coordinates
[pos_S, vel_S] = OE2XCI( satellite_elements_f, GM_PLANET );

% Compute HCI positions of the moon (satellite)
pos = pos + pos_S;
vel = vel + vel_S;
hci = [pos ; vel];

end