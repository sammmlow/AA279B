function E = T2E(v, e)
% Converts true anomaly to eccentric anomaly
%
% Inputs:
%           v - true anomaly (radians)
%           e - eccentricity
% 
% Outputs:
%           E - eccentric anomaly (radians)

v = wrapTo2Pi(v);

% Check for valid eccentricity
if (e >= 0 && e <= 1)

    % Compute eccentric anomaly
    cosE = (e + cos(v)) / (1 + e * cos(v));
    sinE = sqrt(1 - e^2) * sin(v) / (1 + e * cos(v));
    E = atan2(sinE, cosE);
    E = wrapTo2Pi(E);

else
    E = v;
    disp('Invalid eccentricity for T2E')
end