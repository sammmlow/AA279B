function [depArrC3, depSwingbyVminus, swingbyArrVplus, depSwingbyC3] = porkchopPlots(depDates,...
    arrDates, swingbyDates, depPlanetID, swingbyPlanetID, arrPlanetID)

    % Constants:
    muSun = 1.3271244004193938e11;

    % Initialize output arrays:
    nD               = length(depDates);
    nA               = length(arrDates);
    nS               = length(swingbyDates);
    depArrC3         = zeros(nD, nA);
    depSwingbyVminus = zeros(nD, nS);
    swingbyArrVplus  = zeros(nS, nA);
    depSwingbyC3     = zeros(nD, nS);
    
    %%%%%%%%%%%%% DEPARTURE - ARRIVAL: C3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = 1:nD

        rDep = OE2HCI(depPlanetID, depDates(i));

        for j = 1:nA
    
            rArr = OE2HCI(arrPlanetID, arrDates(j));
            tof = (arrDates(j) - depDates(i)) * (24*3600); % seconds
    
            [vdep,~,~] = AA279lambert_curtis(muSun, rDep(1:3), rArr(1:3), 'pro', 0, tof);
    
            % C3:
            v_dep_e_heli  = vdep - rDep(4:6);
            v_dep_e_eci   = rotx(23.45) * v_dep_e_heli;
            depArrC3(i,j) = norm(v_dep_e_eci)^2;
    
        end

    end

    %%%%%%%%%%%%% DEPARTURE - SWINGBY: C3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = 1:nD

        rDep = OE2HCI(depPlanetID, depDates(i));

        for j = 1:nS
    
            rSwing = OE2HCI(swingbyPlanetID, swingbyDates(j));
            tof    = (swingbyDates(j) - depDates(i)) * (24*3600); % seconds
    
            [vdep,~,~] = AA279lambert_curtis(muSun, rDep(1:3), rSwing(1:3), 'pro', 0, tof);
    
            % C3:
            v_dep_e_heli  = vdep - rDep(4:6);
            v_dep_e_eci   = rotx(23.45) * v_dep_e_heli;
            depSwingbyC3(i,j) = norm(v_dep_e_eci)^2;
    
        end

    end

    %%%%%%%%%%%% DEPARTURE - SWINGBY : Swingby Relative Arrival Speed %%%%%

    for i = 1:nD

        rDep = OE2HCI(depPlanetID, depDates(i));

        for j = 1:nS

            rSwing = OE2HCI(swingbyPlanetID, swingbyDates(j));
            tof  = (swingbyDates(j) - depDates(i)) * (24*3600);

            [~, vSwingArr, ~] = AA279lambert_curtis(muSun, rDep(1:3), rSwing(1:3), 'pro', 0, tof);

            % Swingby relative velocity
            depSwingbyVminus(i,j) = norm(vSwingArr - rSwing(4:6));

        end

    end


    %%%%%%%% SWINGBY - ARRIVAL: Swingby Relative Departure Speed %%%%%%%%%%

    for i = 1:nS

        rSwing = OE2HCI(swingbyPlanetID, swingbyDates(i));

        for j = 1:nA

            rArr = OE2HCI(arrPlanetID, arrDates(j));
            tof  = (arrDates(j) - swingbyDates(i)) * (24*3600);

            [vSwingDep,~, ~] = AA279lambert_curtis(muSun, rSwing(1:3), rArr(1:3), 'pro', 0, tof);

            % Swingby relative velocity
            swingbyArrVplus(i,j) = norm(vSwingDep - rSwing(4:6));

        end

    end


end