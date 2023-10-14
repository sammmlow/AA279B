altitude = 1000; % initial altitude in km
musun = 132712440041.93938; % [km^3/sec^2]
muearth = 398600.435436; % [km^3/sec^2] 
muvenus = 324858.592; % [km^3/sec^2]
a_earth = 149597887.5; % km
a_venus = 0.723 * a_earth; % km
ac = 0.5*(a_earth+a_venus);
n_venus = deg2rad( 1.854341616704701E-05 );
n_earth = deg2rad( 1.140795601851852e-05 );

% Compute mean argument of latitudes of planets from JPL Horizons DE441
L_venus = deg2rad( mod( 101.468 + 147.615, 360 ) );
L_earth = deg2rad( mod( 124.262 + 28.383, 360 ) );

% Based on relative period, compute desired period of Hohmann transfer.
a_satellite_helios = a_earth + 6378.140 + altitude;
T_flight = pi * sqrt( ac^3 / musun )
T_flight_years = T_flight / (365.25 * 86400)
DV_initial = sqrt( musun / (a_satellite_helios) );
DV_final = sqrt( musun * ( 2/(a_satellite_helios) - (1/ac) ) );
DV_estimate = DV_final - DV_initial

% Compute the angles (radians) and the wait time.
alpha = T_flight * n_venus
beta_i = L_venus - L_earth
beta_f = pi - alpha

% Compute the wait time
T_wait = (beta_f-beta_i) / (n_venus-n_earth)
T_wait_days = T_wait / 86400
