% T87 atmospheric drag study of Titan 
% R is the height above Titan surface

function rho = titan_density(R)

if R < 980.1
    rho_i = 4.28E-10; % kg/m^3
    Hi = 108.4; % km
    h = 980.1;
elseif R < 1008.4
    rho_i = 4.28E-10; % kg/m^3
    Hi = 108.4; % km
    h = 1008.4;
elseif R < 1045.2
    rho_i = 4.28E-10; % kg/m^3
    Hi = 108.4; % km
    h = 1045.2;
elseif R < 1133.1
    rho_i = 4.28E-10; % kg/m^3
    Hi = 108.4; % km
    h = 1133.1;
else
    rho_i = 0; % kg/m^3
    Hi = 1; % km
    h = R;
end
rho = rho_i * exp( ( R - h ) / Hi );

end