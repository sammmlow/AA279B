function PositionVelocity = computeTrajectory( trajectoryCode, dates )

%{
This function calculates the departure and arrival velocities of a
satellite for each leg of an interplanetary trajectory. 

INPUTS: 
    trajectoryCode - cell array of capital letters corresponding to each 
                     flyby body. Example: Earth-Jupiter-Saturn trajectory 
                     would be {'E','J','S'}.
    dates          - datetime object for each of the expected departure 
                     or arrival to each planet

OUTPUTS:
    N x 12 array where N is the number of transfers. i.e. for 'EJS'
    trajectory, there are two transfers Earth-Jupiter and Jupiter-Saturn 
    so N = 2. For each transfer, the position of the departure and arrival
    planet and the heliocentric velocity of the satellite at the departure 
    and arrival planet is provided:
    [depPlanetPos(1x3), arrPlanetPos(1x3), depVel(1x3), arrVel(1x3)]
%} 
    
    % Initialize Output:
    PositionVelocity = zeros(length(trajectoryCode)-1, 12);
    
    % Check if there is valid number of dates for each leg of trajectory.
    if (length(dates) ~= length(trajectoryCode))
        warning('Number of dates must match number of trajectory legs!');
        PositionVelocity = nan(length(trajectoryCode)-1, 12);
        return;
    end

    keys   = {'Me','V','E','M','J','S','U','N','L','T'};

    for i = 1 : length(trajectoryCode) - 1

        % Planet IDs:
        p1 = trajectoryCode{i};
        p2 = trajectoryCode{i+1};

        % Check if the input keys are invalid?
        if (~any(strcmp(keys, p1)))
            warning('Bad trajectory key! Exiting');
            PositionVelocity = nan(length(trajectoryCode)-1, 12);
            return;
        end

        % Positions:
        p1_eph = ephemeris( p1, dates(i)   );
        p2_eph = ephemeris( p2, dates(i+1) );

        % Arrival and departure velocities:
        tof = seconds(dates(i+1) - dates(i));
        [vdep_hci, varr_hci, ~] = lambert( gravParams('Su') , ...
            p1_eph(1:3), p2_eph(1:3), 'pro', tof);
        
        % Add to output:
        PositionVelocity(i,:) = ...
            [p1_eph(1:3)', p2_eph(1:3)', vdep_hci', varr_hci'];

    end
end
