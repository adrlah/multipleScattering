function TAP_blockageSingleScattering
% % Scenarios for testing blockage and back-scattering
scenario.number = 4; % Scenario from the paper Fig. 4
% scenario.number = 5; % Scenario from the paper Fig. 5

%% Print scenario
fprintf(['\n Scenario ' num2str(scenario.number) ' \n \n'])

%% Implemented model: string
models.List = {'3DFresnel'};

%% PO data
PO.available = 1;

%% MoM data
if scenario.number == 4
    MoM.available = 1;
else
    MoM.available = 0;
end

%% Main code
addpath('Fresnel_Integrals')
scenario                     = getScenarioInfo(scenario);
[txSurf0, rxSurf0, scenario] = getAntennas(scenario);
[scattSurf0, scenario]       = getScatteringSurfaces(scenario);
scenario                     = getSweepVariable(scenario);
PO                           = getPOresults(scenario, PO);
MoM                          = getMoMresults(scenario, MoM);
[scenario, models]           = calculateModels(scenario, models, scattSurf0, txSurf0, rxSurf0);
plotResults(scenario, PO, MoM, models);
end

%% Frequency
function scenario = getScenarioInfo(scenario)
scenario.f      = 28e9;
scenario.c0     = 299792458;
scenario.z0     = 377;
scenario.lambda = scenario.c0/scenario.f;
scenario.k0     = 2*pi/scenario.lambda;
end

%% Antennas
function [txSurf0, rxSurf0, scenario] = getAntennas(scenario)
switch scenario.number
    case 4
        yAntennaWallDistance = 0;
        xAntennaSeparation   = 0;
        zAntHeight           = 0;
        TxPhiAntToWall = pi/2;                         % Point along positive y-axis
        RxPhiAntToWall = TxPhiAntToWall;
        
        txSurf0 = isotropic_TAP();
        rxSurf0 = nardaV637_TAP(scenario.f);
        
        txSurf0.geomT{end+1} = [-xAntennaSeparation/2 -yAntennaWallDistance zAntHeight]; % Move x meters along negative y-axis
        txSurf0.geomT{end+1} = [0 0 1 TxPhiAntToWall]; % Phi rotation
        txSurf0.geomT{end+1} = [-1 0 0];               % Move 1 meter along negative y-axis
        
        rxSurf0.geomT{end+1} = [xAntennaSeparation/2 -yAntennaWallDistance zAntHeight]; % Move x meters along negative y-axis
        rxSurf0.geomT{end+1} = [0 0 1 RxPhiAntToWall]; % Phi rotation
        
    case 5
        TxPhiAntToWall = atan2(5*sin(pi/4), 2*5*cos(pi/4));
        RxPhiAntToWall = pi-atan2(5*sin(pi/4), 2*5*cos(pi/4));
        
        txSurf0              = schwarzbeck9170_TAP(scenario.f);
        rxSurf0              = nardaV637_TAP(scenario.f);
        txSurf0.geomT{end+1} = [-2*5*cos(pi/4) 0 0];  % Move 5 meters along negative _global_ x-axis (local positive y-axis after phiAntToWall rotation)
        rxSurf0.geomT{end+1} = [2*5*cos(pi/4) 0 0];   % Move 10 meters along positive _global_ x/y-axis (after phiAntToWall rotation)
        txSurf0.geomT{end+1} = [0 0 1 TxPhiAntToWall];
        rxSurf0.geomT{end+1} = [0 0 1 RxPhiAntToWall];
        
    otherwise
        error('Scenario not defined!')
end

[~, xyzTx0] = reduce_geomT(txSurf0.geomT);
[~, xyzRx0] = reduce_geomT(rxSurf0.geomT);

scenario.xyzTx0 = xyzTx0;
scenario.xyzRx0 = xyzRx0;

end

%% Scattering surfaces
function [scattSurf0, scenario] = getScatteringSurfaces(scenario)
nWidth    = 2;
nHeight   = 2;
nSamples  = [nWidth nHeight];

switch scenario.number
    case 4
        scattSurf0 = smallRectangleForBlockage(nSamples, scenario.f);
        
    case 5
        w          = 1;
        h          = 1;
        scattSurf0 = rectangleInXzPlane(nSamples, w, h, scenario.f);
        if any(scenario.number == 5)
            scattSurf0(2) = rectangleInXzPlane(nSamples, w, h, scenario.f);
        end
end

%% Gamma
switch scenario.number
    case 4
        gamma = -1; % PEC
    case 5
        gamma = -1*ones(2,1); % PEC
end

scenario.gamma = gamma;
end

%% Sweep variable
function scenario = getSweepVariable(scenario)

resolution = getResolution(scenario);
interval   = getInterval(scenario);

nSamples  = abs((interval(2) - interval(1)))/resolution + 1;
sweepVals = linspace(interval(1), interval(2), nSamples);

idxPO = 1:numel(sweepVals);

if scenario.number == 4
    sweepVals = sweepVals*pi/180;
end

scenario.sweepVals = sweepVals;
scenario.idxPO     = idxPO;

end

function resolution = getResolution(scenario)
switch scenario.number
    case 4 % Rx antenna position in the circle
        resolution = 0.1;  % [degree]
    case 5 % Offset along x-axis
        resolution = 0.01; % [meter]
    otherwise
        resolution = [];
end
end

function interval = getInterval(scenario)
switch scenario.number
    case 4 % Rx antenna position in the circle
        interval = [-90 270];                  % [degree]
    case 5 % Offset along x-axis
        interval = [-5*cos(pi/4) 5*cos(pi/4)]; % [meter]
    otherwise
        interval = [];
end
end

%% Models calculation
function [scenario, models] = calculateModels(scenario, models, scattSurf0, txSurf0, rxSurf0)

% Variable preallocation for speed
scenario.interactionList = [];
models.resultsTh = zeros(numel(scenario.sweepVals), 1);
models.resultsPh = zeros(numel(scenario.sweepVals), 1);
gammaFresnel = zeros(numel(scenario.sweepVals),1);
[gammaFresnelSingleReflection] = deal(gammaFresnel);

if scenario.number == 5
    scattSurf{1} = TAP_moveObject(scattSurf0{1}, [0 5*sin(pi/4) 0]);
end
for k = 1:numel(scenario.sweepVals)  % Move/rotate of surface(s) or antenna (s);
    
    % Move/rotate surface(s)
    switch scenario.number
        case 4
            scattSurf{1} = scattSurf0{1};
            R = 1; % 1 meter radius circle
            x_position        = R*cos(scenario.sweepVals(k));
            y_position        = R*sin(scenario.sweepVals(k));       % y0 = -1 is the starting position
            rxSurf0.geomT{4}  = [x_position y_position 0];
            scenario.xyzRx0   = rxSurf0.geomT{4};
            rxAntennaPointing = atan2(y_position, x_position) - pi; % Point to the center of coordinates
            rxSurf0.geomT{5}  = [0 0 1 rxAntennaPointing];
        case 5
            xyOffset          = [scenario.sweepVals(k) 5*sin(pi/4)/2 0];
            startFromSurface  = 2;
            
            for m = startFromSurface:numel(scattSurf0)
                scattSurf{m}  = TAP_moveObject(scattSurf0{m}, xyOffset);
            end
    end
    
    scattSurfOrdered = scattSurf;
    
    % Coordinates for surfaces needed for Fresnel integrals
    for z = 1:numel(scattSurfOrdered)
        iteration.xyz0(z,:) = scattSurfOrdered{z}.centroid;         % Center of surface
        iteration.xyz1(z,:) = scattSurfOrdered{z}.vertices(1,:);    % Corner 1 of surface
        iteration.t1(z,:)   = scattSurfOrdered{z}.vertices(2,:) - scattSurfOrdered{z}.vertices(1,:);
        iteration.t2(z,:)   = scattSurfOrdered{z}.vertices(3,:) - scattSurfOrdered{z}.vertices(1,:);
        iteration.n(z,:)    = scattSurfOrdered{z}.normalVectors(1,:); % Surface normal
        iteration.area(z,:) = scattSurfOrdered{z}.area;
    end
    
    blockage = iteration;
    
    for m = 1:numel(scattSurfOrdered)
        
        if size(scenario.interactionList,1) < numel(scattSurfOrdered)
            scenario.interactionList = [scenario.interactionList; ...
                num2str(m) blanks(numel(scattSurfOrdered)-numel(num2str(numel(scattSurfOrdered))))];
        end
        
        iteration.xyzTx(m,:) = scenario.xyzTx0;
        iteration.xyzRx(m,:) = scenario.xyzRx0;
        [iteration.xyzRP(m,:), iteration.xyzTxImage(m,:), scenario.blockedPath(m)] = getSpecularReflectionPoint(scenario.xyzTx0, scenario.xyzRx0, iteration.xyz1(m,:), iteration.n(m,:));
        iteration.xyzRxImage(m,:) = getImageLocation(scenario.xyzRx0, iteration.xyz1(m,:), iteration.n(m,:));
        
        % Get distance through specular reflection point
        R0reflex        = norm(scenario.xyzTx0-iteration.xyzRP(m,:)) + norm(scenario.xyzRx0-iteration.xyzRP(m,:));
        scenario.R(k,1) = norm(scenario.xyzRx0-iteration.xyzTxImage(m,:));
        
        % Find point used for free space path loss and antenna gain calculation
        % Point on surface closest to specular reflection point
        [iteration.xyzAmp,s,t] = getPointOnSurfaceClosestToRP(iteration.xyzRP(m,:), iteration.xyz1(m,:), iteration.t1(m,:), iteration.t2(m,:));
        
        % Find end points for Fresnel integrals
        % Project the reflection point on (the extensions of) the four
        % rectangle edges and use these for the calculation "excess distance".
        iteration.xyzEdges          = getProjectionsOnEdges(iteration.xyzRP(m,:), iteration.xyz1(m,:), iteration.t1(m,:), iteration.t2(m,:));
        iteration.xyzEdgesMidpoints = getEdgesMidpoint(scattSurf, iteration.t1(m,:), iteration.t2(m,:));
        
        % Get the gain toward the single ray scattering point
        [iteration.eThTx, iteration.ePhTx, iteration.eThRx, iteration.ePhRx] = ...
            getAntennasTowardSingleRayScatteringPoint(iteration.xyzAmp, scenario.xyzTx0, scenario.xyzRx0, txSurf0, rxSurf0);
        
        % LoS component
        scenario.rLOS(k) = max(0,norm(scenario.xyzRx0-scenario.xyzTx0));
        gffLOS(k)        = exp(-1j*scenario.k0*scenario.rLOS(k))/scenario.rLOS(k) * scenario.lambda/(4*pi); %#ok<*AGROW>
        [eThTx(k), ePhTx(k), eThRx(k), ePhRx(k)] = getAntennaTowardAntenna(scenario.xyzTx0, scenario.xyzRx0, txSurf0, rxSurf0);
        
        % Signs for summing of integrals
        iteration.signs = 1 - 2*[s>0 s<1 t>0 t<1];
        
        % Excess path lengths for Fresnel integrals
        iteration.deltaR(:,m)  = max(0,sqrt(sum(abs2(iteration.xyzEdges-scenario.xyzTx0),2)) + sqrt(sum(abs2(iteration.xyzEdges-scenario.xyzRx0),2)) - R0reflex);
        
        % The phase reference is the specular reflection point.
        rMiddle = norm(iteration.xyzRP(m,:)-scenario.xyzTx0) + norm(iteration.xyzRP(m,:)-scenario.xyzRx0);
        scenario.rAmp(k,1) = norm(iteration.xyzAmp-scenario.xyzTx0) + norm(iteration.xyzAmp-scenario.xyzRx0);
        
        % Green's function times "Friis factor"
        iteration.gff    = exp(-1j*scenario.k0*rMiddle)/rMiddle * scenario.lambda/(4*pi);
        
        % Calculate the 3D Fresnel gamma for a given surface
        gammaFresnelSingleReflection(k,m) = ourFresnelModel(scenario, iteration, m);
        
        % Check if any other surface is blocking the Tx-to-Rx path through the given surface
        if numel(scattSurfOrdered) > 1
            
            
            [blockageMatrixFromTx, blockageMatrixToRx] = checkPathBlockage(scenario, iteration, blockage, m);
            
            if any(blockageMatrixFromTx) || any(blockageMatrixToRx)
                blockage.blockedPath(m) = true;
            else
                blockage.blockedPath(m) = false;
            end
            
            if any(blockageMatrixFromTx)
                
                for n = 1:size(blockageMatrixFromTx, 2)
                    
                    if blockageMatrixFromTx(1, n)
                        
                        blockage.xyzTxImage   = getImageLocation(scenario.xyzTx0, blockage.xyz1(n,:), blockage.n(n,:));
                        blockage.xyzRxImage   = getImageLocation(scenario.xyzRx0, blockage.xyz1(m,:), blockage.n(m,:));
                        [blockage.xyzRP(n,:)] = getSpecularReflectionPoint(blockage.xyzTxImage, blockage.xyzRxImage, blockage.xyz1(n,:), blockage.n(n,:));
                        
                        % Get distance through specular reflection point
                        R0block  = norm(blockage.xyzTxImage-blockage.xyzRP(n,:))+norm(blockage.xyzRxImage-blockage.xyzRP(n,:));
                        
                        % Find end points for Fresnel integrals
                        % Project the reflection point on (the extensions of) the four
                        % rectangle edges and use these for the calculation "excess distance".
                        blockage.xyzEdges = getProjectionsOnEdges(blockage.xyzRP(n,:), blockage.xyz1(n,:), blockage.t1(n,:), blockage.t2(n,:));
                        
                        [blockage.xyzAmp,s,t] = getPointOnSurfaceClosestToRP(blockage.xyzRP(n,:), blockage.xyz1(n,:), blockage.t1(n,:), blockage.t2(n,:));
                        
                        % Signs for summing of integrals
                        blockage.signs = 1 - 2*[s>0 s<1 t>0 t<1];
                        
                        % Excess path lengths for Fresnel integrals
                        blockage.deltaR(:,n) = max(0,sqrt(sum(abs2(blockage.xyzEdges-blockage.xyzTxImage),2)) + sqrt(sum(abs2(blockage.xyzEdges-blockage.xyzRxImage),2)) - R0block);
                        
                        gammaFresnelBlockFromTx(k,n) = ourFresnelModel(scenario, blockage, n);
                    end
                end
            else
                gammaFresnelBlockFromTx(k,:) = zeros(1, size(blockageMatrixFromTx, 2));
            end
            
            if any(blockageMatrixToRx)
                
                for n = 1:size(blockageMatrixToRx,2)
                    
                    if blockageMatrixToRx(1, n)
                        
                        blockage.xyzTxImage   = getImageLocation(scenario.xyzTx0, blockage.xyz1(m,:), blockage.n(m,:));
                        blockage.xyzTxImage   = getImageLocation(blockage.xyzTxImage, blockage.xyz1(n,:), blockage.n(n,:));
                        [blockage.xyzRP(n,:)] = getSpecularReflectionPoint(blockage.xyzTxImage, scenario.xyzRx0, blockage.xyz1(n,:), blockage.n(n,:));
                        
                        % Get distance through specular reflection point
                        R0block = norm(blockage.xyzTxImage-blockage.xyzRP(n,:))+norm(scenario.xyzRx0-blockage.xyzRP(n,:));
                        
                        % Find end points for Fresnel integrals
                        % Project the reflection point on (the extensions of) the four
                        % rectangle edges and use these for the calculation "excess distance".
                        blockage.xyzEdges = getProjectionsOnEdges(blockage.xyzRP(n,:), blockage.xyz1(n,:), blockage.t1(n,:), blockage.t2(n,:));
                        
                        % Find point used for free space path loss and antenna gain calculation
                        [blockage.xyzAmp,s,t] = getPointOnSurfaceClosestToRP(blockage.xyzRP(n,:), blockage.xyz1(n,:), blockage.t1(n,:), blockage.t2(n,:));
                        
                        
                        % Signs for summing of integrals
                        blockage.signs = 1 - 2*[s>0 s<1 t>0 t<1];
                        
                        % Excess path lengths for Fresnel integrals
                        blockage.deltaR(:,n) = max(0,sqrt(sum(abs2(blockage.xyzEdges-blockage.xyzTxImage),2)) + sqrt(sum(abs2(blockage.xyzEdges-scenario.xyzRx0),2)) - R0block);
                        
                        gammaFresnelBlockToRx(k,n) = ourFresnelModel(scenario, blockage, n);
                    end
                end
            else
                gammaFresnelBlockToRx(k,:) = zeros(1, size(blockageMatrixToRx, 2));
            end
            
            if blockage.blockedPath(m)
                gammaFresnel(k,m) = scenario.gamma(m)*gammaFresnelSingleReflection(k,m)*...
                    (1 - sum(gammaFresnelBlockFromTx(k, blockageMatrixFromTx)) ...
                    - sum(gammaFresnelBlockToRx(k, blockageMatrixToRx)));
            else
                gammaFresnel(k,m) = scenario.gamma(m)*gammaFresnelSingleReflection(k,m);
            end
        else
            gammaFresnel(k,m) = scenario.gamma(m)*gammaFresnelSingleReflection(k,m);
        end
        
        modelsTmp(k,m) = gammaFresnel(k,m);
        
        modelsTh(k,m) = modelsTmp(k,m)*iteration.eThTx*iteration.eThRx*iteration.gff;
        modelsPh(k,m) = modelsTmp(k,m)*iteration.ePhTx*iteration.ePhRx*iteration.gff;
    end
end

% LoS component
if scenario.number == 4
    LOScomponentTh = eThTx.*eThRx.*gffLOS;
    LOScomponentPh = ePhTx.*ePhRx.*gffLOS;
end

% Sum the contributions of all interactions.
models.resultsTh = sum(modelsTh, 2);
models.resultsPh = sum(modelsPh, 2);

if scenario.number == 4
    models.resultsTh = models.resultsTh + LOScomponentTh.';
    models.resultsPh = models.resultsPh + LOScomponentPh.';
end
end

function PO = getPOresults(scenario, PO)
tempPO    = loadPOdata(scenario.number);
POresults = tempPO(:);

PO.results = POresults;
end

function MoM = getMoMresults(scenario, MoM)

if MoM.available == 1
    tempMoM = loadMoMdata(scenario.number);
    MoMresults = tempMoM(:);
else
    MoMresults = [];
end

MoM.results = MoMresults;
end

function [scenario, PO, MoM, models] = plotResults(scenario, PO, MoM, models)

if scenario.number == 4
    scenario.sweepVals = scenario.sweepVals/pi*180;
end

% Amplitude plot
[models.tikzAmp, PO.tikzAmp, MoM.tikzAmp]       = plotAmplitudResults(scenario, models, PO, MoM);
% Phase plot
[models.tikzPhase, PO.tikzPhase, MoM.tikzPhase] = plotPhaseResults(scenario, models, PO, MoM);

if PO.available && not(any(strcmp(models.List, 'PO'))) % Add PO to the legend before display
    models.List = [models.List; 'PO'];
end
if MoM.available && not(any(strcmp(models.List, 'MoM'))) % Add MoM to the legend before display
    models.List = [models.List; 'MoM'];
end

models.legend = models.List;

if scenario.number == 4
    figure(999);  legend(models.List)
    figure(1000); legend(models.List)
else
    figure(999);  legend(models.List, 'Location', 'northeastoutside')
    figure(1000); legend(models.List, 'Location', 'northeastoutside')
end
end

function [modelsTikzAmp, POtikzAmp, MoMtikzAmp] = plotAmplitudResults(scenario, models, PO, MoM)

figure(999);
clf

results = models.resultsTh;

xValue = scenario.sweepVals;
modelsTikzAmp = zeros(size(results));
index = 1:numel(xValue);

for iii = 1:numel(models.List)
    modelsTikzAmp(:, iii) = 10*log10(abs2(results(index,iii)));
    if scenario.number == 5
        plot(xValue, modelsTikzAmp(:, iii), 'LineWidth', 1); hold on;
    end
end

if PO.available
    POtikzAmp = 10*log10(abs2(PO.results(index,1)));
    if scenario.number == 5
       plot(xValue, POtikzAmp, 'k--', 'LineWidth', 1); hold on;
    end
else
    POtikzAmp = [];
end

if MoM.available
    MoMtikzAmp = 10*log10(abs2(MoM.results(index,1)));
    if scenario.number == 5
        plot(xValue, MoMtikzAmp, 'm--', 'LineWidth', 1); hold on;
    end
else
    MoMtikzAmp = [];
end

if scenario.number == 4
   plotScenario4Amplitude(scenario, xValue, modelsTikzAmp, POtikzAmp, MoMtikzAmp);
    return
end
getProperAxisLimitsAmplitude(scenario)
getProperXaxisLabel(scenario)
ylabel('Magnitude (dB)')
grid on
end

function [modelsTikzPhase, POtikzPhase, MoMtikzPhase] = plotPhaseResults(scenario, models, PO, MoM)

figure(1000);
clf

results = models.resultsTh;

xValue = scenario.sweepVals;
modelsTikzPhase = zeros(size(results));

index = 1:numel(xValue);

for iii = 1:numel(models.List)
    switch scenario.number
        case 4
            phase_offset = MoM.results(index,1);
        otherwise
            phase_offset = 1*ones(length(index));
    end
    
    if  scenario.number == 4
        modelsTikzPhase(:, iii) = unwrap(angle(results(index,iii)./phase_offset(index,1)))*180/pi;
    else
        modelsTikzPhase(:, iii) = wrapTo360(angle(results(index,iii)./phase_offset(index,1))*180/pi);
    end
    if scenario.number == 5        
        plot(xValue, modelsTikzPhase(:, iii), 'LineWidth', 1); hold on;
    end
end

if PO.available
    if  scenario.number == 4
        POtikzPhase = unwrap(angle(PO.results(index,1)./phase_offset(index,1)))*180/pi;
    else
        POtikzPhase = wrapTo360(angle(PO.results(index,1)./phase_offset(index,1))*180/pi);
    end
    if scenario.number == 5
        plot(xValue, POtikzPhase, 'k--', 'LineWidth', 1); hold on;
    end
else
    POtikzPhase = [];
end

if MoM.available
    if  scenario.number == 4
        MoMtikzPhase = unwrap(angle(MoM.results(index,1)./phase_offset(index,1)))*180/pi;
    else
        MoMtikzPhase = wrapTo360(angle(MoM.results(index,1)./phase_offset(index,1))*180/pi);
    end
    if scenario.number == 5 
        plot(xValue, MoMtikzPhase, 'm--', 'LineWidth', 1); hold on;
    end
else
    MoMtikzPhase = [];
end

if scenario.number == 4
   plotScenario4Phase(scenario, xValue, modelsTikzPhase, POtikzPhase, MoMtikzPhase);
    return
end

getProperAxisLimitsPhase(scenario)
getProperXaxisLabel(scenario);
ylabel('Phase (deg)')
grid on
end

function plotScenario4Amplitude(scenario, x, models, PO, MoM)
subplot(1,3,1)
index = find(x >= -80 & x <= -70);
plot(x(index), models(index, :), 'LineWidth', 1); hold on; grid on;
plot(x(index), PO(index, :), 'k--','LineWidth', 1); hold on; grid on;
plot(x(index), MoM(index, :), 'm--', 'LineWidth', 1); hold on; grid on;
getProperAxisLimitsAmplitude(scenario)
xlim([-80 -70])
getProperXaxisLabel(scenario)
ylabel('Magnitude (dB)')
grid on

subplot(1,3,2)
index = find(x >= -5 & x <= 5);
plot(x(index), models(index, :), 'LineWidth', 1); hold on; grid on;
plot(x(index), PO(index, :), 'k--','LineWidth', 1); hold on; grid on;
plot(x(index), MoM(index, :), 'm--', 'LineWidth', 1); hold on; grid on;
getProperAxisLimitsAmplitude(scenario)
xlim([-5 5])
getProperXaxisLabel(scenario)
ylabel('Magnitude (dB)')
grid on

subplot(1,3,3)
index = find(x >= 80 & x <= 90);
plot(x(index), models(index, :), 'LineWidth', 1); hold on; grid on;
plot(x(index), PO(index, :), 'k--','LineWidth', 1); hold on; grid on;
plot(x(index), MoM(index, :), 'm--', 'LineWidth', 1); hold on; grid on;
getProperAxisLimitsAmplitude(scenario)
xlim([80 90])
getProperXaxisLabel(scenario)
ylabel('Magnitude (dB)')
grid on
end

function plotScenario4Phase(scenario, x, models, PO, MoM)
subplot(1,3,1)
index = find(x >= -80 & x <= -70);
plot(x(index), models(index, :), 'LineWidth', 1); hold on; grid on;
plot(x(index), PO(index, :), 'k--','LineWidth', 1); hold on; grid on;
plot(x(index), MoM(index, :), 'm--', 'LineWidth', 1); hold on; grid on;
getProperAxisLimitsPhase(scenario)
xlim([-80 -70])
getProperXaxisLabel(scenario)
ylabel('Phase (deg)')
grid on

subplot(1,3,2)
index = find(x >= -5 & x <= 5);
plot(x(index), models(index, :), 'LineWidth', 1); hold on; grid on;
plot(x(index), PO(index, :), 'k--','LineWidth', 1); hold on; grid on;
plot(x(index), MoM(index, :), 'm--', 'LineWidth', 1); hold on; grid on;
getProperAxisLimitsPhase(scenario)
xlim([-5 5])
getProperXaxisLabel(scenario)
ylabel('Phase (deg)')
grid on

subplot(1,3,3)
index = find(x >= 80 & x <= 90);
plot(x(index), models(index, :), 'LineWidth', 1); hold on; grid on;
plot(x(index), PO(index, :), 'k--','LineWidth', 1); hold on; grid on;
plot(x(index), MoM(index, :), 'm--', 'LineWidth', 1); hold on; grid on;
getProperAxisLimitsPhase(scenario)
xlim([80 90])
getProperXaxisLabel(scenario)
ylabel('Phase (deg)')
grid on
end


function getProperXaxisLabel(scenario)
switch scenario.number
    case 4
        xlabel('\alpha (deg)')
    case 5
        xlabel('Surface 1 centroid x-coordinate (m)')
    otherwise
        xlabel('Undefined label. Please update')
end
end

function getProperAxisLimitsAmplitude(scenario)
switch scenario.number
    case 4
        ylim([-90 -45]);
    case 5
        xlim([-5*cos(pi/4) 5*cos(pi/4)]);        
        ylim([-70 -45]);       
    otherwise
end
end

function getProperAxisLimitsPhase(scenario)
switch scenario.number
    case 4
        xlim([-90 90]);
        ylim([-90 90]);
    case 5
        xlim([-5*cos(pi/4) 5*cos(pi/4)]);        
        ylim([60 240]); 
    otherwise
end
end

function [POresults] = loadPOdata(scenario)
load(['SimulationData\POscenario' num2str(scenario)], 'POresults');
end

function [MoMresults] = loadMoMdata(scenario)
load(['SimulationData\MoMscenario' num2str(scenario)], 'MoMresults');
end

function xyzTxImage  = getImageLocation(xyzTx,xyzOnSurface,n)
% getImageLocation  Image source in plane containing
% 'xyzOnSurface' with normal 'n'.
xyzTxImage = xyzTx + 2*((xyzOnSurface-xyzTx)*n.')*n; % Image of Tx
end

function [xyzRP, xyzTxImage, blockedPath]  = getSpecularReflectionPoint(xyzTx, xyzRx, xyzOnSurface, n)
% getSpecularReflectionPoint  Specular reflection point on plane containing
% 'xyzOnSurface' with normal 'n'.

blockedPath = false;

% First check whether the Tx and Rx are in different sides of the plane
% defined by the scattering surface

xyzTxImage = xyzTx + 2*((xyzOnSurface-xyzTx)*n.')*n; % Image of Tx
u         = (xyzRx-xyzTxImage)/norm(xyzRx-xyzTxImage);
w         = xyzTxImage-xyzOnSurface;
s         = -(n*w.')/(n*u.');
xyzRP     = xyzTxImage + s*u;

% Check if the segment u1 intersects with the scattering surface. If it
% does, it's a reflection situation; otherwise, a blockage situation
u1         = (xyzRx-xyzTxImage);
s1         = -(n*w.')/(n*u1.');

if (s1 < 0 || s1 > 1) && not(any(ismembertol([0 1], s1))) && not(isinf(s1))
    
    blockedPath = true;
    
    u         = (xyzRx-xyzTx)/norm(xyzRx-xyzTx);
    w         = xyzTx-xyzOnSurface;
    s         = -(n*w.')/(n*u.');
    xyzRP     = xyzTx + s*u;
    xyzTxImage = xyzTx;
    
    return
end
end

function  [blockageMatrixFromTx, blockageMatrixToRx] = checkPathBlockage(scenario, iteration, blockage, interaction)

for k = 1:size(interaction, 2)
    for m = 1:size(blockage.xyz0, 1)
        if any(m == interaction(k))
            blockageMatrixFromTx(k,m) = false;
            blockageMatrixToRx(k,m) = false;
            continue;
        else
            %
            if k == 1
                xyzTx = scenario.xyzTx0;
            else
                xyzTx = iteration.xyzRP(k-1,:);
            end
            if k == size(interaction, 2)
                xyzRx = scenario.xyzRx0;
            else
                xyzRx = iteration.xyzRP(k+1,:);
            end
            
            xyzOnSurface = blockage.xyz1(m,:);
            n = blockage.n(m,:);
            
            xyzTxRP = iteration.xyzRP(k,:);
            xyzRxRP = iteration.xyzRP(k,:);
            
            u1         = (xyzRxRP-xyzTx);
            w1         = xyzTx-xyzOnSurface;
            s1         = -(n*w1.')/(n*u1.');
            
            u2         = (xyzRx-xyzTxRP);
            w2         = xyzOnSurface-xyzRx;
            s2         = -(n*w2.')/(n*u2.');
            
            if (s1 > 0 && s1 < 1) && not(any(ismembertol([0 1], s1)))
                blockageMatrixFromTx(k,m) = true;
            else
                blockageMatrixFromTx(k,m) = false;
            end
            if (s2 > 0 && s2 < 1) && not(any(ismembertol([0 1], s2)))
                blockageMatrixToRx(k,m) = true;
            else
                blockageMatrixToRx(k,m) = false;
            end
        end
    end
end
end

function [xyz0,s0,t0] = getPointOnSurfaceClosestToRP(xyzRp,xyz1moved,t1,t2)
u  = t1/norm(t1); % Unit vector along first edge from first to second corner
v  = t2/norm(t2); % Unit vector along second edge from first to third corner
s0 = (xyzRp-xyz1moved)*u.'/norm(t1); % Projection of vector from rectangle center to reflection point on 'u'
t0 = (xyzRp-xyz1moved)*v.'/norm(t2); % Projection of vector from rectangle center to reflection point on 'v'
s  = s0;
if s<0
    s = 0; % Move to first edge; RP outside rectangle
elseif s>1
    s = 1; % Move to edge opposite of first edge; RP outside rectangle
end
t = t0;
if t<0
    t = 0; % Move to second edge
elseif t>1
    t = 1; % Move to edge opposite of second edge; RP outside rectangle
end
xyz0 = xyz1moved + s*t1 + t*t2;
end

function xyzEdges = getProjectionsOnEdges(xyzRp,xyz1moved,t1,t2)
% Points on line including xyzRp and parallel to to t1 and t2
u = t1/norm(t1); % Unit vector along first edge from first to second corner
v = t2/norm(t2); % Unit vector along second edge from first to third corner

% Project the vector xyz1moved-xyzRp on 'u' to get 'pu', and then
% compute xyzRp + pu and xyzRp + pu + t1 which are the points used for
% the Fresnel integral in the first "dimension" (of the surface). And
% correspondingly for the second dimension.
pu = ((xyz1moved-xyzRp)*u(:))*u; % Projection of xyz1moved-xyzRp on u
pv = ((xyz1moved-xyzRp)*v(:))*v; % Projection of xyz1moved-xyzRp on v
xyzEdges = xyzRp + [pu;pu+t1;pv;pv+t2];
end

function xyzEdgesMidpoints = getEdgesMidpoint(scattSurf,t1,t2)
xyzEdgesMidpoints(1, :) = scattSurf{1}.vertices(1,:) + t2/2;
xyzEdgesMidpoints(2, :) = scattSurf{1}.vertices(2,:) + t2/2;
xyzEdgesMidpoints(3, :) = scattSurf{1}.vertices(1,:) + t1/2;
xyzEdgesMidpoints(4, :) = scattSurf{1}.vertices(3,:) + t1/2;
end

function res = smallRectangleForBlockage(nSamples, freq)
% Test case for comparison with MoM

surfaceWidth  = 0.3;
surfaceHeight = 0.5;
surfaceOffset = [0 0 0];
surfScatt     = TAP_CreateRectangle(surfaceWidth,surfaceHeight,nSamples,freq);
surfScatt     = TAP_rotateObject(surfScatt, [1 0 0], pi/2); % Scattering currents in xz-plane
surfScatt     = TAP_moveObject(surfScatt,surfaceOffset);
res           = {surfScatt};

end

function res = rectangleInXzPlane(nSamples, w, h, freq)
% Rectangle centered at origin
surfScatt     = TAP_CreateRectangle(w,h,nSamples,freq);
surfScatt     = TAP_rotateObject(surfScatt, [1 0 0], pi/2); % Scattering currents in xz-plane
res           = {surfScatt};
end

function [eThTx, ePhTx, eThRx, ePhRx] = getAntennasTowardSingleRayScatteringPoint(xyzAmp, xyzTx0, xyzRx0, txSurf0, rxSurf0)
% Antenna gains toward antenna gain reference point.

thTxMiddle = atan2(sqrt(sum(abs2(xyzAmp(:,1:2)-xyzTx0(1:2)),2)),xyzAmp(:,3)-xyzTx0(3));
phTxMiddle = atan2(xyzAmp(2)-xyzTx0(2),xyzAmp(1)-xyzTx0(1));
thRxMiddle = atan2(sqrt(sum(abs2(xyzAmp(:,1:2)-xyzRx0(1:2)),2)),xyzAmp(:,3)-xyzRx0(3));
phRxMiddle = atan2(xyzAmp(2)-xyzRx0(2),xyzAmp(1)-xyzRx0(1));

rHat_x  = sin(thTxMiddle)*cos(phTxMiddle);
rHat_y  = sin(thTxMiddle)*sin(phTxMiddle);
rHat_z  = cos(thTxMiddle);
R       = reduce_geomT(txSurf0.geomT);
rotRhat = R'*[rHat_x.'; rHat_y.'; rHat_z.'];
[phi,theta] = cart2sph(rotRhat(1,:),rotRhat(2,:),rotRhat(3,:));
theta=pi/2-theta;

[eThTx, ePhTx] = txSurf0.function_handle(theta, phi);

rHat_x  = sin(thRxMiddle)*cos(phRxMiddle);
rHat_y  = sin(thRxMiddle)*sin(phRxMiddle);
rHat_z  = cos(thRxMiddle);
R       = reduce_geomT(rxSurf0.geomT);
rotRhat = R'*[rHat_x.'; rHat_y.'; rHat_z.'];
[phi,theta] = cart2sph(rotRhat(1,:),rotRhat(2,:),rotRhat(3,:));
theta=pi/2-theta;
[eThRx, ePhRx] = rxSurf0.function_handle(theta, phi);

end

function [eThTx, ePhTx, eThRx, ePhRx] = getAntennaTowardAntenna(xyzTx0, xyzRx0, txSurf0, rxSurf0)

thTxToRx   = atan2(sqrt(sum(abs2(xyzRx0(:,1:2)-xyzTx0(1:2)),2)),xyzRx0(:,3)-xyzTx0(3));
phTxToRx   = atan2(xyzRx0(2)-xyzTx0(2),xyzRx0(1)-xyzTx0(1));

thRxToTx   = atan2(sqrt(sum(abs2(xyzTx0(:,1:2)-xyzRx0(1:2)),2)),xyzTx0(:,3)-xyzRx0(3));
phRxToTx   = atan2(xyzTx0(2)-xyzRx0(2),xyzTx0(1)-xyzRx0(1));

rHat_x  = sin(thTxToRx)*cos(phTxToRx);
rHat_y  = sin(thTxToRx)*sin(phTxToRx);
rHat_z  = cos(thTxToRx);
R       = reduce_geomT(txSurf0.geomT);
rotRhat = R'*[rHat_x.'; rHat_y.'; rHat_z.'];
[phi,theta] = cart2sph(rotRhat(1,:),rotRhat(2,:),rotRhat(3,:));
theta=pi/2-theta;
[eThTx, ePhTx] = txSurf0.function_handle(theta, phi);

rHat_x  = sin(thRxToTx)*cos(phRxToTx);
rHat_y  = sin(thRxToTx)*sin(phRxToTx);
rHat_z  = cos(thRxToTx);
R       = reduce_geomT(rxSurf0.geomT);
rotRhat = R'*[rHat_x.'; rHat_y.'; rHat_z.'];
[phi,theta] = cart2sph(rotRhat(1,:),rotRhat(2,:),rotRhat(3,:));
theta=pi/2-theta;
[eThRx, ePhRx] = rxSurf0.function_handle(theta, phi);
end

function gammaFresnel = ourFresnelModel(scenario, iteration, m)
% 3D Fresnel model [1]
% [1] A. Lahuerta-Lavieja, M. Johansson, U. Gustavsson, T. A. H. Bressner, and G. A. E. Vandenbosch, 
% "Computationally-efficient millimeter-wave back-scattering models," to be published in IEEE Trans. Antennas Propag., 2020.

lambda = scenario.lambda;
deltaR = iteration.deltaR(:,m);
signs  = iteration.signs;
fresnelInts = fresnelC(sqrt(4/lambda*deltaR)) - 1j*fresnelS(sqrt(4/lambda*deltaR));
gammaFresnel = 0.5*1j*(((fresnelInts(1)*signs(1)+fresnelInts(2)*signs(2))*(fresnelInts(3)*signs(3)+fresnelInts(4)*signs(4))));
end

function result = abs2(x)
result = real(x).^2 + imag(x).^2;
end

function [surface, nVerts] = TAP_CreateRectangle(width,height,nVerts, frequency)
x     = linspace(0,width,nVerts(1));
x     = x - mean(x);
y     = linspace(0,height,nVerts(2));
y     = y - mean(y);
[x,y] = ndgrid(x,y);
surface.vertices = [x(:) y(:) zeros(size(x(:)))];
nVertsTot        = nVerts(1)*nVerts(2);
faces            = (1:(nVertsTot-nVerts(1))).' + [0 1 nVerts(1)+1 nVerts(1)];
faces(nVerts(1):nVerts(1):end,:) = [];
surface.faces    = faces;
surface.tag = 'rectangleSurface';
surface.frequency=frequency;
surface.type = 'rectangular surface';
siz = size(surface.faces);
surface.centroid = reshape(mean(reshape(surface.vertices(surface.faces.',:),[siz(2) siz(1) 3])),[siz(1) 3]);

faces       = surface.faces;
xV = reshape(surface.vertices(faces,1),siz);
yV = reshape(surface.vertices(faces,2),siz);
zV = reshape(surface.vertices(faces,3),siz);

if size(xV,2)==4
    v1 = [xV(:,3) - xV(:,2) yV(:,3) - yV(:,2) zV(:,3) - zV(:,2)];
    v2 = [xV(:,4) - xV(:,3) yV(:,4) - yV(:,3) zV(:,4) - zV(:,3)];
elseif size(xV,2)==3
    v1 = [xV(:,1) - xV(:,1) yV(:,2) - yV(:,1) zV(:,2) - zV(:,1)];
    v2 = [xV(:,3) - xV(:,2) yV(:,3) - yV(:,2) zV(:,3) - zV(:,2)];
end
xp = [v1(:,2).*v2(:,3)-v2(:,2).*v1(:,3) v1(:,3).*v2(:,1)-v2(:,3).*v1(:,1) v1(:,1).*v2(:,2)-v2(:,1).*v1(:,2)];
nV = xp./sqrt(sum(xp.^2,2));
surface.normalVectors = nV;

crossProd           = zeros([size(surface.faces,1),3,size(surface.faces,2)-2]);
numVerts              = size(surface.faces,2);
for k = 1:numVerts
    k2 = mod(k,numVerts)+1;
    crossProd(:,:,k) = cross([xV(:,k ) yV(:,k ) zV(:,k )],...
        [xV(:,k2) yV(:,k2) zV(:,k2)],2);
end
crossProdSum = sum(crossProd,3);
surface.area  = abs(0.5*sum(surface.normalVectors.*crossProdSum,2));

end

function obj = TAP_moveObject(obj, dist)
switch obj.type
    case 'antenna'
        obj.geomT = {dist obj.geomT{:}};
    otherwise
        obj.vertices = bsxfun(@plus, obj.vertices, dist);
        obj.centroid = bsxfun(@plus, obj.centroid, dist);
end
end

function obj = TAP_rotateObject(obj, rotAxis, alpha)
if alpha==0
    rotMatrixR =eye(3);
else
    rotAxis=rotAxis/sqrt(rotAxis*rotAxis');
    A=[0 -rotAxis(3) rotAxis(2);rotAxis(3) 0 -rotAxis(1);-rotAxis(2) rotAxis(1) 0];
    rotMatrixR = eye(3) + A^2 + A*sin(alpha) - A^2*cos(alpha);
    rotMatrixR = rotMatrixR';
end
switch obj.type
    case 'antenna'
        obj.geomT = {[rotAxis alpha] obj.geomT{:}};
    otherwise
        obj.vertices = obj.vertices*rotMatrixR;
        if isfield(obj, 'centroid')
            obj.centroid = obj.centroid*rotMatrixR;
        end
        if isfield(obj, 'normalVectors')
            obj.normalVectors = obj.normalVectors*rotMatrixR;
        end
end
end

function [R,T] = reduce_geomT(geomT)
affineT = eye(4);
for m = 1:numel(geomT)
    thisT = eye(4);
    if numel(geomT{m}) == 3
        thisT(1:3,4) = geomT{m};
    elseif numel(geomT{m}) == 4
        if geomT{m}(4)==0
            thisT(1:3,1:3)=eye(3);
        else
            n=geomT{m}(1:3);
            n=n/sqrt(n*n');
            A=[0 -n(3) n(2);n(3) 0 -n(1);-n(2) n(1) 0];
            thisT(1:3,1:3) = eye(3) + A^2 + A*sin(geomT{m}(4)) - A^2*cos(geomT{m}(4));
        end
    end
    affineT = affineT*thisT;
end
R = affineT(1:3,1:3);
T = affineT(1:3,4)';
end

function source = nardaV637_TAP(f,varargin)
% nardaV637_TAP
% Leuven Rx horn: Narda V637
% https://nardamiteq.com/product-spec/waveguidehornantennas_Standard_Gain_Horns_2.60_to_40_GHz.pdf
f_scaling               = f/28e9;
para.f                  = f;                  % Frequency (Hz)
% Values and estimates from data sheet and data sheet picture
para.a                  = 0.007112/f_scaling; % Feeding wave guide width (m); WR28
para.b                  = 0.003556/f_scaling; % Feeding wave guide height (m); "-
para.width              = 0.02292/f_scaling;  % Aperture width (m); with wall thickness removed (2x2 mm)
para.height             = 0.01683/f_scaling;  % Aperture height (m); with wall thickness removed (2x2 mm)
para.length             = 0.033724/f_scaling; % Length (along the wave guide feed) of the flaring section
para.scaleToDirectivity = [];
[~,~,para]              = horn_pattern(pi/2,0,para);
fhandle = eval([' @(theta,phi)(',func2str(@horn_pattern),'(theta,phi,para));']);
source.type             = 'antenna';
source.function_handle  = fhandle;
source.geomT            = {};
source.tag              = 'nardav637';
source.geomT{1}      = [0 0 0];
source.geomT{2}      = [1 0 0 0];
if nargin>1 && not(isempty(varargin{1}))
    source.geomT{3} = [0 0 1 varargin{1}];
end
end

function source = schwarzbeck9170_TAP(f,varargin)
% Schwarzbeck9170
% Leuven Tx horn: Schwarzbeck 9170
% http://schwarzbeck.de/Datenblatt/k9170.pdf
f_scaling               = f/28e9;
para.f                  = f;     % Frequency (Hz)
% Values and estimates from data sheet and data sheet picture
para.a                  = 0.010668/f_scaling; % Feeding wave guide width (m); assumed WR42
para.b                  = 0.004318/f_scaling; % Feeding wave guide height (m); "-
para.width              = 0.059/f_scaling;    % Aperture width (m); wall thickness removed
para.height             = 0.044/f_scaling;    % Aperture height (m); wall thickness removed
para.length             = 0.062/f_scaling;    % Length (along the wave guide feed) of the flaring section
para.scaleToDirectivity = [];
[~,~,para]              = horn_pattern(pi/2,0,para);
fhandle = eval([' @(theta,phi)(',func2str(@horn_pattern),'(theta,phi,para));']);
source.type             = 'antenna';
source.function_handle  = fhandle;
source.geomT            = {};
source.tag              = 'schwarzbeck9170';
source.geomT{1} = [0 0 0];
source.geomT{2} = [1 0 0 0];
if nargin>1 && not(isempty(varargin{1}))
    source.geomT{3} = [0 0 1 varargin{1}];
end
end

function [E_theta,E_phi,modelParameters] = horn_pattern(theta,phi,varargin)

if not(isempty(varargin{1}))
    f = varargin{1}.f;
    b = varargin{1}.b;
    a = varargin{1}.a;
    height = varargin{1}.height;
    width = varargin{1}.width;
    length = varargin{1}.length;
    scaleToDirectivity = varargin{1}.scaleToDirectivity;
end

modelParameters.f = f;
modelParameters.b = b;
modelParameters.a = a;
modelParameters.height = height;
modelParameters.width = width;
modelParameters.length = length;
modelParameters.scaleToDirectivity = scaleToDirectivity;

if isempty(scaleToDirectivity)
    tmp = modelParameters;
    tmp.scaleToDirectivity = 1;
    th = linspace(0,pi,361);
    ph = linspace(-pi,pi,721);
    dO = dOmega_TAP(th,ph,[0 pi -pi pi]);
    [th,ph]     = ndgrid(th,ph);
    [E_th,E_ph] = horn_pattern(th,ph,tmp);
    average     = dO*(abs(E_th(:)).^2+abs(E_ph(:)).^2)/4/pi;
    scaleToDirectivity = sqrt(1./average);
    modelParameters.scaleToDirectivity = scaleToDirectivity;
end

[E_theta,E_phi] = alg_horn(theta,phi,f,a,b,width,height,length,scaleToDirectivity);
end

function [E_theta,E_phi] = alg_horn(theta,phi,f,a,b,width,height,length,scaleToDirectivity)
c_0 = 299792458;
k  = 2*pi*f/c_0;
ky = k*sin(theta).*sin(phi);
kz = k*cos(theta);

apexE  = height*sqrt((length/(height-b))^2+.25);
sqrt_factorE = 1./sqrt(pi*apexE*k);

t1      = sqrt_factorE*(-.5*k*height-kz*apexE);
t2      = sqrt_factorE*(.5*k*height-kz*apexE);
fresnel = conj((fresnelC(t2)+1j*fresnelS(t2))-(fresnelC(t1)+1j*fresnelS(t1)));
I2      = sqrt(pi*apexE./k).*exp(.5j*kz.^2.*apexE./k).*fresnel;

apexH  = width*sqrt((length/(width-a))^2+.25);
sqrt_factorH = 1./sqrt(pi*apexH*k);

kyprim = ky+pi/width;
t1prim = sqrt_factorH*(-.5*k*width-kyprim*apexH);
t2prim = sqrt_factorH*(.5*k*width-kyprim*apexH);
fresnelprim = conj((fresnelC(t2prim)+1j*fresnelS(t2prim))-(fresnelC(t1prim)+1j*fresnelS(t1prim)));

kybis  = ky-pi/width;
t1bis = sqrt_factorH*(-.5*k*width-kybis*apexH);
t2bis = sqrt_factorH*(.5*k*width-kybis*apexH);
fresnelbis = conj((fresnelC(t2bis)+1j*fresnelS(t2bis))-(fresnelC(t1bis)+1j*fresnelS(t1bis)));
I1 = .5*sqrt(pi*apexH./k).*(...
    exp(.5j*kyprim.^2.*apexH./k).*fresnelprim + ...
    exp(.5j*kybis.^2.*apexH./k).*fresnelbis);

Amp     = I1.*I2*scaleToDirectivity;
E_theta = -Amp.*(cos(phi) + sin(theta));
E_phi   = Amp.*cos(theta).*sin(phi);
end

function dO=dOmega_TAP(theta,phi,varargin)
limits = varargin{1};
w=(2*mod(floor(theta/pi),2)-1).*cos(theta)+2*floor(theta/pi);
w_av=conv(w,[.5 .5]);
w_av(1)=(2*mod(floor(limits(1)/pi),2)-1).*cos(limits(1))+2*floor(limits(1)/pi);
w_av(end)=(2*mod(floor(limits(2)/pi),2)-1).*cos(limits(2))+2*floor(limits(2)/pi);
dw=abs(diff(w_av));
phi_av=conv(phi,[.5 .5]);
phi_av(1)=limits(3);
phi_av(end)=limits(4);
dphi=abs(diff(phi_av));
[dw,dphi]=ndgrid(dw,dphi);
dO=(dw(:).*dphi(:))';
end

function source = isotropic_TAP(varargin)
parameters.pol = [1 0];
fhandle = eval([' @(theta,phi)(',func2str(@isotropic_pattern),'(theta,phi,parameters));']);
source.type             = 'antenna';
source.function_handle  = fhandle;
source.geomT            = {};
source.tag              = 'isotropic';
source.geomT{1}         = [0 0 0];
source.geomT{2}         = [1 0 0 0];
if nargin>0 && not(isempty(varargin{1}))
    source.geomT{3}     = [0 0 1 varargin{1}];
end
end

function [E_theta, E_phi, modelParameters] = isotropic_pattern(theta, phi, varargin)
pol  = [1 0];
modelParameters.pol = pol;
E_theta=pol(1)*ones(size(theta+phi));
E_phi  =pol(2)*ones(size(theta+phi));
end
