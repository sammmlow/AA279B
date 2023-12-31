%% Generate Porkchop Plot Arrays:

% test script to see if I can generate an unpowered interplanetary
% trajectory for Earth-Jupiter-Saturn

% Constant
jdOffset = 1721058.5;

earthDepInit = datetime(2035,6,15,0,0,0);
satArrInit   = datetime(2040,7,5,0,0,0);
depDates     = juliandate(earthDepInit): 2 : juliandate(earthDepInit + 365*3);
arrDates     = juliandate(satArrInit)  : 5 : juliandate(satArrInit + 365*6);
swingbyDates = depDates(1) + 500       : 5 : arrDates(1) - 100;

[depArr, depSwingby, swingbyArr, depSwingbyC3] = porkchopPlots(depDates,...
    arrDates, swingbyDates, 3, 5, 6);

%% Porkchop plots

depDatesDateStr     = datestr(depDates - jdOffset);
arrDatesDateStr     = datestr(arrDates - jdOffset);
swingbyDatesDateStr = datestr(swingbyDates - jdOffset);

figure()
contourLevels = 100:10:200;
contourf(depDates - jdOffset, arrDates - jdOffset,...
    depArr', contourLevels,'showtext','on','LineWidth',2, 'FaceAlpha', 0.5)
% contourf(depDatesDateStr' , arrDatesDateStr',...
%     depArr', contourLevels,'showtext','on','LineWidth',2, 'FaceAlpha', 0.5, 'LabelSpacing', 250)
title('C3 for Earth-Saturn')
% xticklabels({depDatesDateStr})
% yticklabels({arrDatesDateStr})
datetick('x', 'mmm-yyyy');
datetick('y', 'mmm-yyyy');
% xlabel(['Earth Launch Date [days after ' char(earthDepInit) ']'])
% ylabel(['Saturn Arrival Date [days after ' char(earthDepInit) ']'])
xlabel('Earth Launch Date')
ylabel('Saturn Arrival Date')
fontsize(gcf,18, "points")
grid on;

figure()
% contourLevels = round(min(depSwingby,[],'all')) : 2 : round(max(depSwingby,[],'all'));
contourLevels = 0:40;
contourf(swingbyDates - jdOffset, depDates - jdOffset,...
    depSwingby,contourLevels,'showtext','on','LineWidth', 2, 'FaceAlpha', 0.5, 'LabelSpacing', 250);
title('Jupiter Arrival Speed for Earth-Jupiter')
ylabel('Earth Launch Date')
xlabel('Jupiter Swingby Date')
datetick('x', 'mmm-yyyy');
datetick('y', 'mmm-yyyy');
fontsize(gcf,18, "points")
grid on;

figure()
contourLevels = round(min(swingbyArr,[],'all')) : 2 : round(max(swingbyArr,[],'all'));
contourf(swingbyDates - jdOffset, arrDates - jdOffset,...
    swingbyArr', contourLevels,'showtext','on','LineWidth',2, 'FaceAlpha', 0.5, 'LabelSpacing', 250)
title('Jupiter Departure Speed for Jupiter-Saturn')
ylabel('Saturn Arrival Date')
xlabel('Jupiter Swingby Date')
datetick('x', 'mmm-yyyy');
datetick('y', 'mmm-yyyy');
fontsize(gcf,18, "points")
grid on;

figure()
% contourLevels = round(min(depSwingbyC3,[],'all')) : 2 : round(max(depSwingbyC3,[],'all'));
contourLevels = 50:10:300;
contourf(depDates - jdOffset, swingbyDates - jdOffset,...
    depSwingbyC3',contourLevels,'showtext','on','LineWidth',2, 'FaceAlpha', 0.5, 'LabelSpacing', 250)
title('C3 for Earth-Jupiter')
xlabel('Earth Launch Date')
ylabel('Jupiter Swingby Date')
datetick('x', 'mmm-yyyy');
datetick('y', 'mmm-yyyy');
fontsize(gcf,18, "points")
grid on;
%% Create new EJS porkchop plot :

jdOffset = 1721058.5;
muSun    = 1.3271244004193938e11;
radius_J = 69911;

% Create empty array for new EJS porkchop plot (Earth Departure, Saturn
% Arrival, [Earth-Jupiter C3, Jupiter flyby date])
EJS_C3        = zeros(length(depDates), length(arrDates), 2);
EJS_C3(:,:,1) = -999; % Set all date combos as impossible until we find a flyby date that works
EJS_C3(:,:,2) = 0;

% Choose a tolerance for range of C3 values to test above the min C3 value
C3tol = 200;

% Choose the Earth departure - Saturn arrival that minimizes C3
[minValue, minIdx] = min(depArr, [], 'all');
[departureIdx, arrivalIdx] = find(depArr < C3tol);

for j = 1:length(departureIdx)
    % Using the Earth departure - Saturn arrival date, sweep to find the
    % Jupiter flyby that satisfies that combo
    depIdx = departureIdx(j);
    arrIdx = arrivalIdx(j);

    % fprintf("Earth-Saturn C3 for this iteration is: %.4g \n", depArr(depIdx, arrIdx))

    flybyDateIdx = findFlybyDates(departureIdx(j), arrivalIdx(j), depSwingby, swingbyArr, 0.01);

    if ~isempty(flybyDateIdx)

        % If there are flyby dates, test that the flyby dates actually do produce an unpowered
        % flyby:
        flybyC3 = [];
        
        for i = 1:length(flybyDateIdx)

            muJ = 126686511;
            eJD = depDates(depIdx);
            sJD = arrDates(arrIdx);
            jJD = swingbyDates(flybyDateIdx(i));

            if eJD < jJD && jJD < sJD
            
                tofEJ = (jJD - eJD) * (24 * 3600);                    
                rE = OE2HCI(3, eJD);
                rJ = OE2HCI(5, jJD);
    
                [vdep,varr,~] = AA279lambert_curtis(muSun, rE(1:3), rJ(1:3), 'pro', 0, tofEJ);

                % C3:
                v_dep_e_heli     = vdep - rE(4:6);
                v_dep_e_eci      = rotx(23.45) * v_dep_e_heli;
                flybyC3(end + 1) = norm(v_dep_e_eci)^2;
                
    
                % tofJS = (sJD - jJD) * (24 * 3600);
                % rS = OE2HCI(6, sJD);
                % [vdep_hci, varr_hci, vdep_pl,varr_pl, vdep_norm, varr_norm, delta, e, rp]...
                %     = PlanetSwingby(muJ, rE(1:3), rJ(1:3), rS(1:3), tofEJ, tofJS, rJ(4:6));


                % 
                % fprintf(['Earth Departure Date: ' char(datestr(eJD - jdOffset)) '\n'])
                % fprintf(['Jupiter Swingby Date: ' char(datestr(jJD - jdOffset)) '\n'])
                % fprintf(['Saturn Arrival Date: ' char(datestr(sJD - jdOffset)) '\n'])
                % fprintf("Planetocentric arrival and departure speed: %.4g and %.4g km/s \n", varr_norm, vdep_norm);
                % fprintf("Radius of closest approach is: %.4g \n", rp);
                % fprintf("Required Earth-Saturn C3: %.4g km^2/s^2 \n", depArr(depIdx,arrIdx))
                % fprintf("Required Earth-Jupiter C3: %.4g km^2/s^2 \n", flybyC3(end))
                % fprintf("\n")

            else
                continue
            end
        
        end

        if rp > radius_J

            [minVal, minIdx] = min(flybyC3);
            EJS_C3(depIdx, arrIdx, 1) = minVal;
            EJS_C3(depIdx, arrIdx, 2) = swingbyDates(flybyDateIdx(minIdx));
        else
            continue
        end

    else 
        continue;
        
    end
end

%% Earth-Jupiter-Saturn Porkchop Plot

figure()
% contourLevels = linspace(0,max(EJS_C3(:,:,1), [], 'all'), 20);
contourLevels = 50:10:200;
contourf(depDates - jdOffset, arrDates - jdOffset,...
    EJS_C3(:,:,1)', contourLevels, 'showtext','on','LineWidth', 1.5, 'FaceAlpha', 0.5)
title('C3 for Earth-Jupiter-Saturn')
xlabel('Earth Launch Date')
ylabel('Saturn Arrival Date')
datetick('x', 'mmm-yyyy');
datetick('y', 'mmm-yyyy');
grid on;

fontsize(gcf,18,"points")
%% Find all dates corresponding to C3 values beneath a certain threshold:

C3_threshold = 150;
logicalMatrix = (EJS_C3(:,:,1) > 0) & (EJS_C3(:,:,1) < C3_threshold);
[departureIdxAll, arrivalIdxAll] = find(logicalMatrix);
% [departureIdxAll, arrivalIdxAll] = find(EJS_C3(EJS_C3(:,:,1) < C3_threshold & EJS_C3(:,:,1) > 0));

SaturnArrivalAll = zeros(length(departureIdxAll), 7);
DatesAll         = zeros(length(departureIdxAll), 3); %[eJD, jJD, sJD] 

muJ     = 126686511;
AU      = 149597870.7;
muSun   = 1.3271244004193938e11;
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
tstep   = 3600; % Seconds per hour

for k = 1:length(departureIdxAll)

    % Calculate trajectory and add to output array
    
    eJD = depDates(departureIdxAll(k));
    sJD = arrDates(arrivalIdxAll(k));
    jJD = EJS_C3(departureIdxAll(k), arrivalIdxAll(k), 2);

    tofEJ = (jJD - eJD) * (24 * 3600); 
    tofJS = (sJD - jJD) * (24 * 3600);
    rE = OE2HCI(3, eJD);
    rJ = OE2HCI(5, jJD);
    rS = OE2HCI(6, sJD);

    % Earth   -> Jupiter
    tspanEJ        = 0 : tstep : tofEJ; % simulate once an hour
    [vdep_E, ~, ~] = AA279lambert_curtis(muSun, rE(1:3), rJ(1:3), 'pro', 1, tofEJ);
    initialEJ      = [rE(1:3); vdep_E];
    [~,State_EJ]   = ode113(@(t,State) Propogate2Body(State, muSun),...
                     tspanEJ, initialEJ, options);
    % Jupiter -> Saturn
    tspanJS      = 0 : tstep : tofJS; % simulate once an hour
    [vdep_J, ~, ~] = AA279lambert_curtis(muSun, rJ(1:3), rS(1:3), 'pro', 1, tofJS);
    initialJS    = [rJ(1:3); vdep_J];
    [~,State_JS] = ode113(@(t,State) Propogate2Body(State, muSun),...
                     tspanJS, initialJS, options);

    SaturnArrivalAll(k,:) = [sJD, State_JS(end,:)];
    DatesAll(k,:)         = [eJD, jJD, sJD];

    
    % % Sanity Check (suppress this later)
    % [vdep_hci, varr_hci, vdep_pl,varr_pl, vdep_norm, varr_norm, delta, e, rp]...
    % = PlanetSwingby(muJ, rE(1:3), rJ(1:3), rS(1:3), tofEJ, tofJS, rJ(4:6));
    % 
    % fprintf(['Earth Departure Date: ' char(datestr(eJD - jdOffset)) '\n'])
    % fprintf(['Jupiter Swingby Date: ' char(datestr(jJD - jdOffset)) '\n'])
    % fprintf(['Saturn Arrival Date: ' char(datestr(sJD - jdOffset)) '\n'])
    % fprintf("Planetocentric arrival and departure speed: %.4g and %.4g km/s \n", varr_norm, vdep_norm);
    % fprintf("Radius of closest approach is: %.4g \n", rp);
    % fprintf("Required C3: %.4g km^2/s^2 \n", EJS_C3(departureIdxAll(k),arrivalIdxAll(k),1))


end
%% Find the minimum C3 Value and Corresponding julian dates

[minC3, minC3Idx] = min(EJS_C3(EJS_C3 > 0), [], 'all');
[departureIdx, arrivalIdx] = find(EJS_C3(:,:,1) == minC3);

if length(departureIdx) > 1
    departureIdx = departureIdx(1);
    arrivalIdx   = arrivalIdx(1);
end

muJ = 126686511;
eJD = depDates(departureIdx);
sJD = arrDates(arrivalIdx);
jJD = EJS_C3(departureIdx, arrivalIdx, 2);

tofEJ = (jJD - eJD) * (24 * 3600); 
tofJS = (sJD - jJD) * (24 * 3600);
rE = OE2HCI(3, eJD);
rJ = OE2HCI(5, jJD);
rS = OE2HCI(6, sJD);
[vdep_hci, varr_hci, vdep_pl,varr_pl, vdep_norm, varr_norm, delta, e, rp]...
    = PlanetSwingby(muJ, rE(1:3), rJ(1:3), rS(1:3), tofEJ, tofJS, rJ(4:6));

fprintf(['Earth Departure Date: ' char(datestr(eJD - jdOffset)) '\n'])
fprintf(['Jupiter Swingby Date: ' char(datestr(jJD - jdOffset)) '\n'])
fprintf(['Saturn Arrival Date: ' char(datestr(sJD - jdOffset)) '\n'])
fprintf("Planetocentric arrival and departure speed: %.4g and %.4g km/s \n", varr_norm, vdep_norm);
fprintf("Radius of closest approach is: %.4g \n", rp);
fprintf("Required Earth-Jupiter C3: %.4g km^2/s^2 \n", EJS_C3(departureIdx,arrivalIdx,1))
fprintf("Required Earth-Saturn C3 for same date: %.4g km^2/s^2 \n", depArr(departureIdx, arrivalIdx));

%% Plot trajectory:

AU      = 149597870.7;
muSun   = 1.3271244004193938e11;
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
tstep   = 3600; % Seconds per hour

% Earth   -> Jupiter
tspanEJ        = 0 : tstep : tofEJ; % simulate once an hour
[vdep_E, ~, ~] = AA279lambert_curtis(muSun, rE(1:3), rJ(1:3), 'pro', 1, tofEJ);
initialEJ      = [rE(1:3); vdep_E];
[~,State_EJ]   = ode113(@(t,State) Propogate2Body(State, muSun),...
                     tspanEJ, initialEJ, options);
% Jupiter -> Saturn
tspanJS      = 0 : tstep : tofJS; % simulate once an hour
[vdep_J, ~, ~] = AA279lambert_curtis(muSun, rJ(1:3), rS(1:3), 'pro', 1, tofJS);
initialJS    = [rJ(1:3); vdep_J];
[~,State_JS] = ode113(@(t,State) Propogate2Body(State, muSun),...
                     tspanJS, initialJS, options);

figure()
% Orbital Radii
centers = zeros(3,2);
radii   = [1; 5.2; 9.5];
viscircles(centers,radii,'color','black','LineWidth',1.5);
hold on;
% Planet Positions
scatter(rE(2)/AU, rE(1)/AU, 100, 'filled')
text(rE(2)/AU -1, rE(1)/AU +1, string(datestr(eJD - jdOffset)))
scatter(rJ(2)/AU, rJ(1)/AU, 100, 'filled')
text(rJ(2)/AU -1, rJ(1)/AU +1, string(datestr(jJD - jdOffset)))
scatter(rS(2)/AU, rS(1)/AU, 100, 'filled')
text(rS(2)/AU -1, rS(1)/AU +1, string(datestr(sJD - jdOffset)))
% Trajectory
plot(State_EJ(:,2)/AU,State_EJ(:,1)/AU,'r-','LineWidth',2)
plot(State_JS(:,2)/AU,State_JS(:,1)/AU,'r-','LineWidth',2)
xlabel('Y [AU]')
ylabel('X [AU]')
title('Trajectory from Earth to Saturn with Unpowered Jupiter Flyby')
legend('Earth','Jupiter','Saturn')
axis equal;
grid on;
set(gca, 'XDir', 'reverse')
hold off;

%% OUTPUT FOR SAM!!!!

% Output array that is: [saturn_arrival_epoch, rx, ry, rz, vx, vy, vz]
SaturnArrivalBest = [sJD, State_JS(end,:)];








