function [pos, vel] = OE2XCI(OE, MU)
% Take orbital elements and compute position and velocity in ECI frame
%
% Inputs: 
%   oe - 6x1 orbit element vectors consisting of
%        a - semi-major axis (km)
%        e - eccentricity
%        i - inclination (radians)
%        O - right ascension of ascending node (radians)
%        w - argument of periapsis (radians)
%        va - true anomaly (radians)
%   MU - gravitational parameter of planet (km^3/s^2)
% 
% Outputs:
%   pos - 3x1 XCI position state vector [km]
%   vel - 3x1 XCI velocity state vector [km/s]

eps = 1E-5; % Tolerance for integer rounding for e and i

% Check for valid inputs
if size(OE, 1) == 6 && size(OE, 2) == 1 && MU > 0
    
    % Get individual orbit elements
    a = OE(1);
    e = OE(2);
    i = wrapToPi(OE(3));
    O = wrapTo2Pi(OE(4));
    w = wrapTo2Pi(OE(5));
    va = wrapTo2Pi(OE(6));

    % Test for valid inputs
    if (e >= 0 && e <= 1 && a >= 0)

        % Test for edge cases
        if (e < eps)
            % Circular orbit
            if (i < eps) || (abs(i - pi) < eps)
                % Circular equatorial orbit
                w = 0;
                O = 0;
                va = O + w + va;
            else
                % Circular inclined orbit
                w = 0;
                va = w + va;
            end
        else
            if (i < eps) || (abs(i - pi) < eps)
                % Elliptical equatorial orbit
                O = 0;
                w = O + w; 
            end
        end

        % Test for parabolic orbit
        if (abs(e - 1) < eps)
            disp('Parabolic orbit: please input rp in place of a.')
            p = 2 * a;
        else
            p = a * (1 - e^2);
        end
        % Compute radius and velocity of orbit
        r = p / (1 + e * cos(va));
        v = sqrt(MU / p);

        % Find position and velocity in perifocal coordinates
        rpqw = r * [cos(va); sin(va); 0];
        vpqw = v * [-sin(va); e + cos(va); 0];

        % Define rotation matrix from perifocal to ECI (passive rotations)
        % R = Rz(-rasc)Rx(-inc)Rz(-argp)
        R1 = [cos(-O) sin(-O) 0;...
             -sin(-O) cos(-O) 0;...
              0       0       1];
        R2 = [1  0        0;...
              0  cos(-i) sin(-i);...
              0 -sin(-i) cos(-i)];
        R3 = [cos(-w) sin(-w) 0;...
             -sin(-w) cos(-w) 0;...
              0        0        1];
        R = R1 * R2 * R3;
        
        % Find position and velocity in XCI coordinates
        pos = R * rpqw;
        vel = R * vpqw;

    else
        pos = nan(3,1);
        vel = nan(3,1);
        disp('Invalid orbit elements in OE2XCI');
    end

else
    pos = nan(3,1);
    vel = nan(3,1);
    disp('Invalid inputs to OE2XCI');
end