function PositionVelocity = computeTrajectory(trajectoryCode, dates)

%{ 
This function calculates the departure and arrival velocities of a
satellite for each leg of an interplanetary trajectory. 

INPUTS: 
    trajectoryCode - character array of capital letters corresponding to each flyby
                     planet. example: Earth-Jupiter-Saturn trajectory would
                     be 'EJS'
    dates          - datetime object for each the expected departure or arrival
                     to each planet

OUTPUTS:
    N x 12 array where N is the number of transfers. i.e. for 'EJS'
    trajectory, there are two transfers Earth-Jupiter and Jupiter-Saturn so
    N = 2. For each transfer, the position of the departure and arrival
    planet and the heliocentric velocity of the satellite at the departure 
    and arrival planet is provided: [depPlanetPos(1x3), arrPlanetPos(1x3), depVel(1x3), arrVel(1x3)]

%} 

    % dictionary to map planets to planetID:
    keys   = {'Me','V','E','M','J','S','U','N'};
    values = [1, 2, 3, 4, 5, 6, 7, 8];
    planetID = dictionary(keys,values);

    % Initialize Output:
    PositionVelocity = zeros(length(trajectoryCode)-1, 12);

    for i = 1:length(trajectoryCode) - 1

        % Planet IDs:
        p1 = planetID({trajectoryCode(i)});
        p2 = planetID({trajectoryCode(i+1)});

        % Positions:
        p1_pos = OE2HCI(p1,juliandate(dates(i)));
        p2_pos = OE2HCI(p2,juliandate(dates(i+1)));

        % Velocities:
        tof = seconds(dates(i+1) - dates(i));
        [vdep_hci, varr_hci, ~] = AA279lambert_curtis(gravParams({'Su'}),...
            p1_pos(1:3), p2_pos(1:3), 'pro', 1, tof);
        
        % Add to output:
        PositionVelocity(i,:) = [p1_pos(1:3)', p2_pos(1:3)', vdep_hci', varr_hci'];


    end

end
