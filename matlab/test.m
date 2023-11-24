%% Generate Porkchop Plot Arrays:

% Constant
jdOffset = 1721058.5;

earthDepInit = datetime(2037,9,15,0,0,0);
satArrInit   = datetime(2040,7,5,0,0,0);
depDates     = juliandate(earthDepInit): 2 : juliandate(earthDepInit + 120);
arrDates     = juliandate(satArrInit)  : 10 : juliandate(satArrInit + 365*6);
swingbyDates = depDates(1) + 500       : 10 : arrDates(1) - 100;

[depArr, depSwingby, swingbyArr] = porkchopPlots(depDates,...
    arrDates, swingbyDates, 3, 5, 6);

%% Porkchop plots

depDatesDateStr     = datestr(depDates - jdOffset);
arrDatesDateStr     = datestr(arrDates - jdOffset);
swingbyDatesDateStr = datestr(swingbyDates - jdOffset);


figure()
contourLevels = 100:10:200;
contourf(depDates - juliandate(earthDepInit),arrDates - juliandate(earthDepInit),...
    depArr', contourLevels,'showtext','on','LineWidth',2, 'FaceAlpha', 0.5)
% contourf(depDatesDateStr' , arrDatesDateStr',...
%     depArr', contourLevels,'showtext','on','LineWidth',2, 'FaceAlpha', 0.5, 'LabelSpacing', 250)
title('C3 for Earth-Saturn')
% xticklabels({depDatesDateStr})
% yticklabels({arrDatesDateStr})
xlabel(['Launch Date [days after ' char(earthDepInit) ']'])
ylabel(['Arrival Date [days after ' char(earthDepInit) ']'])
grid on;

figure()
contourLevels = round(min(depSwingby,[],'all')) : 2 : round(max(depSwingby,[],'all'));
contourf(depDates - juliandate(earthDepInit),swingbyDates - juliandate(earthDepInit),...
    depSwingby',contourLevels,'showtext','on','LineWidth', 2, 'FaceAlpha', 0.5, 'LabelSpacing', 250);
title('Vminus for Earth-Jupiter')
xlabel(['Launch Date [days after ' char(earthDepInit) ']'])
ylabel(['Jupiter Swingby [days after ' char(earthDepInit) ']'])
grid on;

figure()
contourLevels = round(min(swingbyArr,[],'all')) : 2 : round(max(swingbyArr,[],'all'));
contourf(arrDates - juliandate(earthDepInit),swingbyDates - juliandate(earthDepInit),...
    swingbyArr, contourLevels,'showtext','on','LineWidth',2, 'FaceAlpha', 0.5, 'LabelSpacing', 250)
title('Vpls for Earth-Saturn')
xlabel(['Arrival Date [days after ' char(earthDepInit) ']'])
ylabel(['Jupiter Swingby [days after ' char(earthDepInit) ']'])
grid on;


%% Finding Jupiter flyby dates:

jdOffset = 1721058.5;

[minValue, minIdx] = min(depSwingby, [], 'all');
[departureIdx, arrivalIdx] = find(depSwingby == minValue);

flybyDates = findFlybyDates(departureIdx, arrivalIdx, depSwingby, swingbyArr, 0.05);

% Now let's test that the flyby dates actually do produce an unpowered
% flyby:

for i = 1:length(flybyDates)

    muJ = 126686511;
    eJD = depDates(departureIdx);
    sJD = arrDates(arrivalIdx);
    jJD = swingbyDates(flybyDates(i));
    
    tofEJ = (jJD - eJD) * (24 * 3600); 
    tofJS = (sJD - jJD) * (24 * 3600);
    
    rE = OE2HCI(3, eJD);
    rJ = OE2HCI(5, jJD);
    rS = OE2HCI(6, sJD);
    
    [vdep_hci, varr_hci, vdep_pl,varr_pl, vdep_norm, varr_norm, delta, e, rp]...
        = PlanetSwingby(muJ,rE(1:3), rJ(1:3), rS(1:3), tofEJ, tofJS, rJ(4:6));
    
    fprintf(['Earth Departure Date: ' char(datestr(eJD - jdOffset)) '\n'])
    fprintf(['Jupiter Swingby Date: ' char(datestr(jJD - jdOffset)) '\n'])
    fprintf(['Saturn Arrival Date: ' char(datestr(sJD - jdOffset)) '\n'])
    fprintf("Planetocentric arrival and departure speed: %.4g and %.4g km/s \n", varr_norm, vdep_norm);

end

%% Plot trajectory:












