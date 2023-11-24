function [vdep_hci, varr_hci, vdep_pl,...
    varr_pl, vdep_norm, varr_norm, delta, e, rp] = PlanetSwingby(mu, r1, r2, r3, tof1, tof2, v_planet)

%{ 
This function calculates the asymptotic velocities when approaching and
departing a "swingby" planet both in heliocentric and planetocentric
coordinates. It additionally provides the turning angle and hyperbolic
eccentricity of the swingby trajectory.

INPUTS: 
    mu       - gravitational parameter of the swingby planet (used for
               eccentricity calculation
    r1       - position vector in heliocentric coord of departing planet
    r2       - position vector in heliocentric coord of the swingby planet
    r3       - position vector in heliocentric coord of destination planet
    tof1     - time-of-flight (seconds) of departing to swingby
    tof2     - time-of-flight (seconds) of swingby to destination
    v_planet - heliocentric velocity of swingby planet

OUTPUTS:
    varr_hci  - heliocentric asymptotic velocity arriving at the swingby planet
    vdep_hci  - heliocentric asymptotic velocity departing swingby planet
    varr_pl   - planetocentric asymptotic velocity arriving at the swingby planet
    vdep_pl   - planetocentric asymptotic velocity departing swingby planet
    varr_norm - planetocentric asymptotic speed arriving at the swingby planet
    vdep_norm - planetocentric asymptotic speed departing swingby planet
    delta     - turning angle [rad]
    e         - hyperbolic eccentricity

%} 

 % Heliocentric asymptotic velocities
 muSun = 1.3271244004193938e11;
 [~, varr_hci, ~] = AA279lambert_curtis(muSun, r1, r2, 'pro', 1, tof1);
 [vdep_hci, ~, ~] = AA279lambert_curtis(muSun, r2, r3, 'pro', 1, tof2);

 % Planetocentric asymptotic velocities
 varr_pl = varr_hci - v_planet;
 vdep_pl = vdep_hci - v_planet;
 
 % Planetocentric speeds
 varr_norm = norm(varr_pl);
 vdep_norm = norm(vdep_pl);
 
 % Turning Angle and eccentricity:
 delta = acos(dot(varr_pl,vdep_pl)/(varr_norm*vdep_norm));
 e     = 1/(sin(delta/2));
 rp    = ((e - 1)* mu)/(varr_norm * vdep_norm);


end