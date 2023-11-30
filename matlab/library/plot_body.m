function plot_body( body, radius, coord)
% Inputs:
%   - body   --> char array e.g. 'Jupiter'
%   - radius --> positive scalar float radius in km
%   - coord  --> length 3 vector for coordinates of planet

addpath("../graphics");
hold on;

% Select the image file based on the chosen body
image_file = '';
switch body
    case {'Me'}      % Mercury
        image_file = 'bodyplot_mercury.jpg';
    case {'V'}       % Venus
        image_file = 'bodyplot_venus.jpg';
    case {'E'}       % Earth
        image_file = 'bodyplot_earth.jpg';
    case {'L'}       % Moon
        image_file = 'bodyplot_moon.jpg';
    case {'M'}       % Mars
        image_file = 'bodyplot_mars.jpg';
    case {'J'}       % Jupiter
        image_file = 'bodyplot_jupiter.jpg';
    case {'S'}       % Saturn
        image_file = 'bodyplot_saturn.jpg';
    case {'T'}       % Titan
        image_file = 'bodyplot_venus.jpg';
    otherwise
        warning('Input body does not have an image file! Exiting');
        return;
end

% Create a 3D ellipsoid
[X,Y,Z] = sphere;
X = X * radius + coord(1);
Y = Y * radius + coord(2);
Z = Z * radius + coord(3);
globe = surf(X, Y, -Z, 'FaceColor', 'none', 'EdgeColor', [0 0 0]);

% Texture-map the globe
cdata = imread(image_file);

% Set image as color data (cdata) property, and set face color to
% indicate a texturemap, which Matlab expects to be in cdata.
% Turn off the mesh edges.
set(globe, 'FaceColor', 'texturemap', 'CData', cdata, ...
    'FaceAlpha', 1, 'EdgeColor', 'none');
end