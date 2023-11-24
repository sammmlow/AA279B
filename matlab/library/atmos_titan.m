%  FILE DESCRIPTION:                                                    
%                                                                       
%  Atmospheric density model of Titan using the T87 drag study.
%  Returns the atmospheric density (kg/m^3) given an altitude input (km)
%  Valid only for altitudes up to 1133.1 km.                     
%                                                                       
%  Written by Samuel Y. W. Low.                                         
%  First created 21-Mar-2022                          
%  Last modified 25-Mar-2022

function density = atmos_titan(R)

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
density = rho_i * exp( ( R - h ) / Hi );

end