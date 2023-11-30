%% Function computes a full trajectory from Saturn to Titan, given some
% initial conditions. This function should be run for each gridcell of
% a contour plot, that takes as input a time-of-flight (arrival at Titan),
% an initial epoch, position, and velocity (HCI) at Saturn departure,
% and produces the trajectory. From the trajectory, we can find the
% arrival excess velocity.

close all; clear all; clc;
set(groot, 'defaultLineLineWidth', 1);
addpath("library");

GM_TITAN = 8.9685e3;
GM_SATURN = 3.7941e7;

RAD_TITAN = 2575;
RAD_SATURN = 58232;

%% Function inputs (for use later) all in HCI

% Titan target in the form of orbital elements
epoch = datetime(2045, 8, 28, 0, 0, 0);
tflight = 86400.0 * 2.5;
x0 = 1e9 * [ -0.499580881317805 ; ...
             -1.410423254627437 ; ...
              0.044349388468382 ; ...
              0.000000005976360 ; ...
             -0.000000003709886 ; ...
             -0.000000000054102 ];

%% Construct parameters and transforms between HCI and Synodic frames

eph_saturn = ephemeris( 'S', epoch );
eph_titan = ephemeris( 'T', epoch );
titan_from_saturn = eph_titan - eph_saturn;

% HCI coordinates of Titan relative to Saturn
rT = titan_from_saturn(1:3);
vT = titan_from_saturn(4:6);

% Columns of the direction cosine matrices transforming SYN2HCI
xhat_hci = rT / norm( rT );
vhat_hci = vT / norm( vT );
zhat_hci = cross( xhat_hci, vhat_hci );
yhat_hci = cross( zhat_hci, xhat_hci );

% Form the 3x3 direction cosine matrix between the HCI and Synodic Frame
syn2hci = [ xhat_hci yhat_hci zhat_hci ];
hci2syn = transpose( syn2hci );

% Form the 6x6 HCI to Synodic conversion (apply to position + velocity)
oH = cross(rT, vT) / dot(rT, rT); % Omega expressed in HCI coordinates
omega_hci2syn = [ 0 -oH(3) oH(2) ; oH(3) 0 -oH(1) ; -oH(2) oH(1) 0 ];
rotation_hci2syn = (-1) * hci2syn * omega_hci2syn;
hci2syn6 = [ hci2syn zeros(3,3) ; rotation_hci2syn hci2syn ];

% Form the 6x6 Synodic to HCI conversion (apply to position + velocity)
oS = [0.0 ; 0.0 ; norm(oH)]; % Omega expressed in Synodic coordinates
omega_syn2hci = [ 0 -oS(3) oS(2) ; oS(3) 0 -oS(1) ; -oS(2) oS(1) 0 ];
rotation_syn2hci = syn2hci * omega_syn2hci;
syn2hci6 = [ syn2hci zeros(3,3) ; rotation_syn2hci syn2hci ];

% Scalar (+) distances of Saturn (R1) & Titan (R2) from barycenter [km]
R12 = norm( eph_titan(1:3) - eph_saturn(1:3) );
R1 = -(GM_TITAN  / (GM_TITAN + GM_SATURN)) * R12; 
R2 = (GM_SATURN / (GM_TITAN + GM_SATURN)) * R12; 

% Synodic angular velocity
ws = norm(oH);

% Convert the initial state to synodic frame positions and velocities.
% This involves changing the origin to Saturn, and performing a rotation.
x0_syn = hci2syn6 * (x0 - eph_saturn) + [R1;0;0;0;0;0];
x0_hci_test = syn2hci6 * (x0_syn - [R1;0;0;0;0;0]) + eph_saturn;

%% Now, lets define a conversion from HCI to Titan-Centered Inertial
% tilt = Saturn's inclination + Titan's inclination + Titan's tilt
tilt = 2.486 + 26.73 + 0.34854; 
hci2tci = [1.0         0.0        0.0 ;
           0.0  cosd(tilt) sind(tilt) ;
           0.0 -sind(tilt) cosd(tilt) ];
hci2tci6 = [ hci2tci zeros(3,3) ; zeros(3,3) hci2tci ];

%% Generate the plot object in the Synodic reference frame
figure(1);
scatter3( [R1,R2], [0,0], [0,0], 8, 'black', 'filled');
hold on; grid on; axis equal;
plot_body( 'S', RAD_SATURN, [R1 0 0] );
hold on; grid on;   
plot_body( 'T', RAD_TITAN,  [R2 0 0] );

%% Stage 1: CR3BP shooting method - go for near Titan Center

iter = 1;
maxIters = 20;
maxIters2 = 3 * maxIters;
tstep = 600; % [sec]
error = Inf;
x0_iter = x0_syn;
tolerance = 100.0;
perturbation = 0.001; % [km/s]
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
colormaps = turbo(maxIters2);

target = [ (R2 + RAD_TITAN) ; 0 ; 0];

while error > tolerance

    % For the very first iteration, compute the Jacobian. The Jacobian
    % quantifies the sensitivity of the final distance and speed with 
    % respect to Titan, given some perturbed velocity components in XYZ.
    jacobian = zeros(3,3);
    
    % For each dimension of the velocity vector,
    for idx = 1 : 3

        % The component of that dimension is perturbed,
        perturbed_vector = [ 0 ; 0 ; 0 ; 0 ; 0 ; 0 ];
        perturbed_vector(idx+3) = perturbation;
        
        % One step forward,
        x0_forward = x0_iter + perturbed_vector;
        [tf, xf_forward] = ode113( @(t,y) ...
            derivative(t, y, GM_SATURN, GM_TITAN, R1, R2), ...
            [0:tstep:tflight]', x0_forward, options);

        % One step back...
        x0_backward = x0_iter - perturbed_vector;
        [tf, xf_backward] = ode113( @(t,y) ...
            derivative(t, y, GM_SATURN, GM_TITAN, R1, R2), ...
            [0:tstep:tflight]', x0_backward, options);

        % Compute the Jacobian column.
        jacobian(:,idx) = ...
            (xf_forward(end,1:3)' - xf_backward(end,1:3)') / ...
            (2 * perturbation);
        
    end

    % Apply the shooting method update now.
    [tf, xf_iter] = ode113( @(t,y) ...
            derivative(t, y, GM_SATURN, GM_TITAN, R1, R2), ...
            [0:tstep:tflight]', x0_iter, options);
    target_error = xf_iter(end,1:3)' - target;
    error = norm(target_error);

    % Plot the intermediate trajectories.
    plot3(xf_iter(:,1), xf_iter(:,2), xf_iter(:,3), ...
        'Color', [colormaps(iter,:) 0.75]);

    % Now apply the update to the control variables.
    initial_velocity = x0_iter(4:6);
    updated_velocity = initial_velocity - jacobian \ target_error;
    x0_iter(4:6) = updated_velocity; % Overwrite
    iter = iter + 1;

    disp(['Stage 1: Distance norm error = ' num2str(error)]);
    
    % TODO: Output a NaN if shooting method fails to provide a solution
    if (iter > maxIters)
        break;
    end
end

%% Set x0_iter to some mid-course correction, arbitrarily set to 80% mark
mark = 0.8;
N = length(xf_iter(:,1));
N_mark = round(mark * N);
tflight2 = (1 - mark) * tflight;
x0_iter_mark = xf_iter(N_mark,:)';

%% Stage 2: CR3BP shooting method - now, aim for a circular orbit

err_r = Inf; % [km]
err_i = Inf; % [deg]
err_e = Inf; % [nil]

target_r = 5000; 
target_i = 78.5;
target_e = 0.0;

tol_r = 10;   % [km]
tol_i = 1.0;  % [deg]
tol_e = 0.1;  % [nil]

perturbation = perturbation / 10;

while (err_r > tol_r) || (err_e > tol_e) || (err_i > tol_i)

    % For the very first iteration, compute the Jacobian. The Jacobian
    % quantifies the sensitivity of the final distance and speed with 
    % respect to Titan, given some perturbed velocity components in XYZ.
    jacobian = zeros(3,3);
    
    % For each dimension of the velocity vector,
    for idx = 1 : 3

        % The component of that dimension is perturbed,
        perturbed_vector = [ 0 ; 0 ; 0 ; 0 ; 0 ; 0 ];
        perturbed_vector(idx+3) = perturbation;
        
        % One step forward,
        x0_forward = x0_iter_mark + perturbed_vector;
        [tf, xf_forward] = ode113( @(t,y) ...
            derivative(t, y, GM_SATURN, GM_TITAN, R1, R2), ...
            [0:tstep:tflight2]', x0_forward, options);

        % One step back...
        x0_backward = x0_iter_mark - perturbed_vector;
        [tf, xf_backward] = ode113( @(t,y) ...
            derivative(t, y, GM_SATURN, GM_TITAN, R1, R2), ...
            [0:tstep:tflight2]', x0_backward, options);

        % Convert Synodic to Titan-centered inertial coordinates.
        
        xf_forward_hci = syn2hci6 * (xf_forward(end,:)' - ...
            [R1;0;0;0;0;0]) + eph_saturn;
        xf_forward_tci = hci2tci6 * (xf_forward_hci - eph_titan);

        xf_backward_hci = syn2hci6 * (xf_backward(end,:)' - ...
            [R1;0;0;0;0;0]) + eph_saturn;
        xf_backward_tci = hci2tci6 * (xf_backward_hci - eph_titan);

        % Compute Jacobian with respect to distance to Titan
        j_dist_tci = (norm(xf_forward_tci(1:3)) - ...
            norm(xf_backward_tci(1:3))) / (2 * perturbation);
        
        % Compute Jacobian with respect to eccentricity about Titan
        oe_forward = XCI2OE( xf_forward_tci, GM_TITAN );
        oe_backward = XCI2OE( xf_backward_tci, GM_TITAN );
        e_forward = oe_forward(2);
        e_backward = oe_backward(2);
        j_eccentric = (e_forward - e_backward) / (2 * perturbation);

        % Compute Jacobian with respect to orbit inclination about Titan
        i_forward = rad2deg(oe_forward(3));
        i_backward = rad2deg(oe_backward(3));
        j_inclination = (i_forward - i_backward) / (2 * perturbation);

        % Populate the Jacobian
        jacobian(1,idx) = j_dist_tci;
        jacobian(2,idx) = j_inclination;
        jacobian(3,idx) = j_eccentric;
    end

    % Apply the shooting method update now.
    [tf, xf_iter] = ode113( @(t,y) ...
            derivative(t, y, GM_SATURN, GM_TITAN, R1, R2), ...
            [0:tstep:tflight2]', x0_iter_mark, options);

    xf_iter_hci = syn2hci6 * (xf_iter(end,:)' - ...
            [R1;0;0;0;0;0]) + eph_saturn;
    xf_iter_tci = hci2tci6 * (xf_iter_hci - eph_titan);

    xf_iter_dist = norm(xf_iter_tci(1:3));
    xf_iter_oe = XCI2OE( xf_iter_tci, GM_TITAN );
    xf_iter_incli = rad2deg(xf_iter_oe(3));
    xf_iter_eccen = xf_iter_oe(2);

    err_r = xf_iter_dist  - target_r;
    err_i = xf_iter_incli - target_i;
    err_e = xf_iter_eccen - target_e;

    target_error = [ err_r ; err_i ; err_e ];

    disp(['Stage 2: Error = ' num2str(target_error')]);

    % Artificially dampen the effect of the error
    target_error = target_error * 0.011;

    % Plot the intermediate trajectories.
    plot3(xf_iter(:,1), xf_iter(:,2), xf_iter(:,3), ...
        'Color', [colormaps(iter,:) 0.75]);

    % Now apply the update to the control variables.
    initial_velocity = x0_iter_mark(4:6);
    updated_velocity = initial_velocity - jacobian \ target_error;
    x0_iter_mark(4:6) = updated_velocity; % Overwrite
    iter = iter + 1;
    
    % TODO: Output a NaN if shooting method fails to provide a solution
    if (iter > maxIters2)
        break;
    end
end

%% Stage 3: Check the eccentricity is less than one.
xf_iter_end = xf_iter(end,:);
xf_hci = syn2hci6 * (xf_iter_end' - [R1;0;0;0;0;0]) + eph_saturn;
xf_tci = hci2tci6 * (xf_hci - eph_titan);
xf_tci_r = norm(xf_tci(1:3));
xf_tci_v = norm(xf_tci(4:6));
ev = ((xf_tci_v^2 - GM_TITAN / xf_tci_r) * xf_tci(1:3) - ...
    dot( xf_tci(1:3), xf_tci(4:6) ) * xf_tci(4:6)) / GM_TITAN;
    
% Compute eccentricity
e = norm(ev);
disp(['Eccentricity at arrival to Titan = ' num2str(e)]);

%% Generate legends and labels for the plot.
legend_array = ["","",""];
for n = 1:iter-1
    legend_string = "Iteration " + num2str(n);
    legend_array = [legend_array legend_string];
end
legend([legend_array]);
xlabel('Saturn-Titan Synodic Frame X [km]');
ylabel('Saturn-Titan Synodic Frame Y [km]');

%% Derivative computation for ODE113 (synodic frame trajectory)
% Note: it is assumed that R1 < 0, and R2 > 0
function dsdt = derivative(t, y, mu1, mu2, R1, R2)

    % Preliminary values for clarity.
    w = sqrt(((mu1 + mu2) ) / (-R1+R2)^3);
    pos_R1 = [  R1 ; 0.0 ; 0.0 ];
    pos_R2 = [  R2 ; 0.0 ; 0.0 ];
    pos_R13 = y(1:3) - pos_R1; % Position of S/C w.r.t Saturn
    pos_R23 = y(1:3) - pos_R2; % Position of S/C w.r.t Titan
    acc_R13 = mu1 * pos_R13 / ( norm(pos_R13)^3 );
    acc_R23 = mu2 * pos_R23 / ( norm(pos_R23)^3 );

    % Actual computation of derivatives.
    dsdt = [ y(4) ; ...
             y(5) ; ...
             y(6) ; ...
             2*w*y(5) + w*w*y(1) - acc_R13(1) - acc_R23(1); ...
            -2*w*y(4) + w*w*y(2) - acc_R13(2) - acc_R23(2); ...
            -acc_R13(3) - acc_R23(3)];
end


%% Function to update the Jacobian using Broyden's Method
