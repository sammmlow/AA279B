function oe = XCI2OE(XCI, MU)
% Take Planet-X-Centered Inertial vectors and compute orbital elements
%
% Inputs: 
%   XCI - 6x1 XCI state vectors consisting of three position components 
%         and three velocity components (km, km, km, km/s, km/s, km/s)
%   MU - gravitational parameter of planet (km^3/s^2)
%
% Outputs: 
%   oe - 6x1 orbit element state vectors consisting of
%        a - semi-major axis (km)
%        e - eccentricity
%        i - inclination (radians)
%        O - right ascension of ascending node (radians)
%        w - argument of periapsis (radians)
%        va - true anomaly (radians)

% NOTE THAT THIS FUNCTION HAS BEEN MODIFIED TO OUTPUT HYPERBOLIC PARAMS

% Tolerance for integer rounding for e and i
eps = 1E-5;

% Check for valid inputs
if size(XCI, 1) == 6 && size(XCI, 2) == 1 && MU > 0

    % Size output
    oe = zeros(6,1);
        
    % Extract position and velocity
    reci = XCI(1:3);
    veci = XCI(4:6);
    r = norm(reci);
    v = norm(veci);
    % Orbital angular momentum
    hv = cross(reci, veci);
    h = norm(hv);
    nv = cross([0 0 1], hv);
    n = norm(nv);

    % Compute inclination
    i = acos(hv(3) / h); % Returns values in interval [0,pi]
    i = wrapToPi(i);
    if i < eps
        i = 0;
    elseif abs(i - pi) < eps
        i = pi;
    end

    % Eccentricity vector
    ev = ((v^2 - MU / r) * reci - dot(reci, veci) * veci) / MU;
    
    % Compute eccentricity
    e = norm(ev);

    % Check for valid value
    if (e < 1)
    
        % Specific mechanical energy
        sme = v^2 / 2 - MU / r;
        % Assume integer value if within tolerance eps
        if e < eps
            e = 0;
        elseif abs(1 - e) < eps
            e = 1;
        end

        % Compute a (test for parabolic orbit)
        if e ~= 1
            a = -MU / (2 * sme);
        else
            a = inf;
        end

        % Compute RAAN, argument of periapsis, and mean anomaly
        if (e < eps)
            % Circular equatorial (O, w, v undefined)
            if (i == 0) || (i == pi)
                disp('The orbit is circular equatorial. True longitude is returned in place of O, w and v.')
                truel = reci(1) / r;
                if reci(2) < 0
                    truel = 2*pi - truel;    % Quadrant test
                end
                O = truel;
                w = truel;
                va = truel;
            % Circular inclined (w, v undefined)
            else
                disp('The orbit is circular inclined. Argument of latitude is returned in place of w and v.')
                argl = acos(dot(nv, reci) / (n * r));
                if reci(3) < 0
                    argl = 2*pi - argl;  % Quadrant test
                end
                O = acos(nv(1) / norm(nv));
                if nv(2) < 0
                    O = 2*pi - O;  % Quadrant test
                end
                w = argl;
                va = argl;
            end
        else
            % Elliptical equatorial (O, w undefined)
            if (i == 0) || (i == pi)
                disp('The orbit is elliptical equatorial. True longitude of periapsis is returned in place of O and w.')
                truelp = acos(e(1) / e);
                if ev(2) < 0
                    truelp = 2*pi - truelp;  % Quadrant test
                end
                O = truelp;
                w = truelp;
                va = acosd(dot(ev , reci) / (e * r));
                if dot(reci, veci) < 0
                    va = 2*pi - va;    % Quadrant test
                end
            % Elliptical inclined
            else
                O = acos(nv(1) / norm(nv));
                if nv(2) < 0
                    O = 2*pi - O;  % Quadrant test
                end
                w = acos(dot(nv, ev) / (n * e));
                if ev(3) < 0
                    w = 2*pi - w;  % Quadrant test
                end
                va = acos(dot(ev, reci) / (e * r));
                if dot(reci, veci) < 0
                    va = 2*pi - va;    % Quadrant test
                end
            end
        end 

        % Assign outputs
        oe(1) = a;
        oe(2) = e;
        oe(3) = i;
        oe(4) = O;
        oe(5) = w;
        oe(6) = va;
    else
        oe(1) = r; % This is technically wrong, but just hack it
        oe(2) = e;
        oe(3) = i;
        oe(4) = nan;
        oe(5) = nan;
        oe(6) = nan;
        % disp('Warning! Hyperbolic eccentricity detected!');
    end

else
    oe = nan(6,1);
    disp('Invalid inputs to XCI2OE');
end