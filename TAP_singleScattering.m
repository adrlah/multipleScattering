function TAP_singleScattering
%% Scenarios for testing the single-scattering models
% scenario.number = 2; % Scenario from the paper Fig. 2
% scenario.number = 4; % Scenario from the paper Fig. 4
% scenario.number = 5; % Scenario from the paper Fig. 5
scenario.number = 6; % Scenario from the paper Fig. 6
% scenario.number = 7; % Scenario from the paper Fig. 7

%% Print scenario
fprintf(['\n Scenario ' num2str(scenario.number) ' \n \n'])

%% modelslList: string
% Options: '3DFresnel', 'erf', 'M-METIS', 'M-mmMAGIC', 'M-ITUFresnel', 'M-ITUEmpirical', 'METISRCS'
models.List = {'3DFresnel'; 'erf'; 'M-METIS'; 'M-mmMAGIC';  'M-ITUFresnel'; 'M-ITUEmpirical'};
% models.List = {'3DFresnel'; 'erf'; 'M-METIS'; 'METISRCS'; 'M-mmMAGIC';  'M-ITUFresnel'; 'M-ITUEmpirical';}; % METIS RCS only used in scenario/Figure 7
% models.List = {'3DFresnel'};

%% PO data
PO.available = 1;

%% MoM data
if scenario.number == 2
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
getErrorMetric(scenario, PO, MoM, models);
plotResults(scenario, PO, MoM, models);
end

%% Frequency
function scenario = getScenarioInfo(scenario)
switch scenario.number
    case {4}
        scenario.f = 60e9;
    case {5}
        scenario.f = 39e9;
    otherwise
        scenario.f = 28e9;
end
scenario.c0     = 299792458;
scenario.z0     = 377;
scenario.lambda = scenario.c0/scenario.f;
scenario.k0     = 2*pi/scenario.lambda;
end

%% Antennas
function [txSurf0, rxSurf0, scenario] = getAntennas(scenario)

switch scenario.number
    case {2, 4, 5, 6, 7}
        switch scenario.number
            case 2
                yAntennaWallDistance = 1;
                xAntennaSeparation   = 1;
                zAntHeight           = 0;
            case 4
                yAntennaWallDistance = 1;
                xAntennaSeparation   = 0;
                zAntHeight           = 0;
            case 5
                yAntennaWallDistance = 1;
                xAntennaSeparation   = 2;
                zAntHeight           = 0;
            case 6
                yAntennaWallDistance = 1.3918;
                xAntennaSeparation   = 6.7164;
                zAntHeight           = 1.6;
            otherwise
                yAntennaWallDistance = 0;
                xAntennaSeparation   = 0;
                zAntHeight           = 0;
        end
        switch scenario.number
            case 2
                TxPhiAntToWall       = atan2(yAntennaWallDistance, 0.5*xAntennaSeparation);
                RxPhiAntToWall       = pi;
            case {4, 5}
                TxPhiAntToWall       = atan2(yAntennaWallDistance, 0.5*xAntennaSeparation);
                RxPhiAntToWall       = pi-TxPhiAntToWall;
            case 6
                TxPhiAntToWall       = 1.8437/180*pi;  % Pointing error
                RxPhiAntToWall       = pi;
            case 7
                TxPhiAntToWall       = pi/2;          
                RxPhiAntToWall       = TxPhiAntToWall;
            otherwise
                TxPhiAntToWall       = atan2(yAntennaWallDistance, 0.5*xAntennaSeparation);
                RxPhiAntToWall       = TxPhiAntToWall;
        end
        
        switch scenario.number
            case 2
                txSurf0      = isotropic_TAP();                 
                rxSurf0      = nardaV637_TAP(scenario.f);      
            otherwise
                txSurf0      = schwarzbeck9170_TAP(scenario.f); 
                rxSurf0      = nardaV637_TAP(scenario.f);      
        end
        
        txSurf0.geomT{end+1} = [-xAntennaSeparation/2 -yAntennaWallDistance zAntHeight];
        txSurf0.geomT{end+1} = [0 0 1 TxPhiAntToWall];                                  
        
        if scenario.number == 7
            txSurf0.geomT{end+1} = [-2 0 0];                                            
        end
        
        rxSurf0.geomT{end+1} = [xAntennaSeparation/2 -yAntennaWallDistance zAntHeight]; 
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
    case {2, 7}
        scattSurf0 = rectangleForMoMorDistance(nSamples, scenario.f);
    case {4, 5}
        scattSurf0 = rectangleForRotation(nSamples, scenario.f);
    case 6
        scattSurf0 = rectangleLeuven(nSamples, scenario.f);
    otherwise
end

%% Gamma
switch scenario.number
    otherwise
        gamma = -1; % PEC
end

scenario.gamma = gamma;

end

%% Sweep variable
function scenario = getSweepVariable(scenario)

resolution = getResolution(scenario);
interval   = getInterval(scenario);

if scenario.number == 7
    nSamples = 1e3;
    sweepVals = logspace(interval(1), interval(2), nSamples);
else
    nSamples  = abs((interval(2) - interval(1)))/resolution + 1;
    sweepVals = linspace(interval(1), interval(2), nSamples);
end

idxPO = 1:numel(sweepVals);

if any(scenario.number == [2, 4, 5])
    sweepVals = sweepVals*pi/180;
end

scenario.sweepVals = sweepVals;
scenario.idxPO     = idxPO;

end

function resolution = getResolution(scenario)
switch scenario.number
    case {2}                 % Rx antenna pointing angle
        resolution = 1;      % [degree]
    case {4, 5}              % Rotation in the yz-plane
        resolution = 0.5;    % [degree]
    case {6}                 % Offset along x-axis
        resolution = 0.002;  % [meter]
    otherwise
        resolution = [];
end
end

function interval = getInterval(scenario)
switch scenario.number
    case {2}                    % Rx antenna pointing angle
        interval = [-90 180];   % [degree]
    case {4, 5}                 % Rotation in the yz-plane
        interval = [0 180];     % [meter]
    case 6                      % Offset along x-axis
        interval = [-1 1];      % [meter]
    case {7}                    % Offset along y-axis
        interval = [-1 3];      % [10^(interval) meter]
    otherwise
        interval = [];
end
end

%% Models calculation
function [scenario, models] = calculateModels(scenario, models, scattSurf0, txSurf0, rxSurf0)

% Variable preallocation for speed
scenario.interactionList = [];
modelsTmp = zeros(numel(models.List),numel(scenario.sweepVals));
[modelsTh, modelsPh] = deal(modelsTmp);
models.resultsTh = zeros(numel(scenario.sweepVals), numel(models.List));
models.resultsPh = zeros(numel(scenario.sweepVals), numel(models.List));
gammaFresnel = zeros(numel(scenario.sweepVals),1);
[gammaErf, gammaMetis, gammaMMMagic, gammaMetisRCS, gammaITUFres, gammaITUEmp] = deal(gammaFresnel);
[gammaFresnelSingleReflection, gammaErfSingleReflection, gammaMetisSingleReflection, ...
    gammaMMMagicSingleReflection, gammaMetisRCSSingleReflection, gammaITUFresSingleReflection, gammaITUEmpSingleReflection] = deal(gammaFresnel);

for k = 1:numel(scenario.sweepVals)  % Move/rotate of surface(s) or antenna(s);
    
    switch scenario.number
        case {2}
            scattSurf{1} = scattSurf0{1};
            rxSurf0.geomT{4} = [0 0 1 -(pi+scenario.sweepVals(k))];
        case {4, 5}
            xyzC = scattSurf0{1}(1).centroid;
            scattSurf{1} = moveObject(scattSurf0{1}(1),[-xyzC(1) 0 0]);
            scattSurf{1} = rotateObject(scattSurf{1},[0 1 0], scenario.sweepVals(k));
            scattSurf{1} = moveObject(scattSurf{1},[xyzC(1) 0 0]);
        case 6
            xyOffset  = [scenario.sweepVals(k) 0 0];
            scattSurf{1} = moveObject(scattSurf0{1}, xyOffset);
        case {7}
            xyOffset  = [0 scenario.sweepVals(k) 0];
            clear scattSurf
            scattSurf{1} = moveObject(scattSurf0{1},xyOffset);
        otherwise
    end
    
    scattSurfOrdered = scattSurf;
    
    % Coordinates for surfaces needed for Fresnel integrals
    for z = 1:numel(scattSurfOrdered)
        iteration.xyz0(z,:) = scattSurfOrdered{z}.centroid;         
        iteration.xyz1(z,:) = scattSurfOrdered{z}.vertices(1,:);    
        iteration.t1(z,:)   = scattSurfOrdered{z}.vertices(2,:) - scattSurfOrdered{z}.vertices(1,:);
        iteration.t2(z,:)   = scattSurfOrdered{z}.vertices(3,:) - scattSurfOrdered{z}.vertices(1,:);
        iteration.n(z,:)    = scattSurfOrdered{z}.normalVectors(1,:); 
        iteration.area(z,:) = scattSurfOrdered{z}.area;
    end
    
    iteration.xyzTx(1,:) = scenario.xyzTx0;
    iteration.xyzRx(1,:) = scenario.xyzRx0;
    [iteration.xyzRP(1,:), iteration.xyzTxImage(1,:)] = getSpecularReflectionPoint(scenario.xyzTx0, scenario.xyzRx0, iteration.xyz1(1,:), iteration.n(1,:));
    iteration.xyzRxImage(1,:) = getImageLocation(scenario.xyzRx0, iteration.xyz1(1,:), iteration.n(1,:));
    
    % Get distance through specular reflection point
    R0reflex        = norm(scenario.xyzTx0-iteration.xyzRP(1,:)) + norm(scenario.xyzRx0-iteration.xyzRP(1,:));
    scenario.R(k,1) = norm(scenario.xyzRx0-iteration.xyzTxImage(1,:));
    
    % Find antenna gain reference point used for free space path loss and antenna gain calculation
    [iteration.xyzAmp,s,t] = getPointOnSurfaceClosestToRP(iteration.xyzRP(1,:), iteration.xyz1(1,:), iteration.t1(1,:), iteration.t2(1,:));
    
    % Find end points for Fresnel integrals   
    iteration.xyzEdges          = getProjectionsOnEdges(iteration.xyzRP(1,:), iteration.xyz1(1,:), iteration.t1(1,:), iteration.t2(1,:));
    iteration.xyzEdgesMidpoints = getEdgesMidpoint(scattSurf, iteration.t1(1,:), iteration.t2(1,:));
    
    % Get the gain toward the antenna gain reference point
    [iteration.eThTx, iteration.ePhTx, iteration.eThRx, iteration.ePhRx] = ...
        getAntennasTowardSingleRayScatteringPoint(iteration.xyzAmp, scenario.xyzTx0, scenario.xyzRx0, txSurf0, rxSurf0);
    
    if scenario.number == 2 % Get gain for Tx-to-Rx LoS component
        [~, ~, g.eThRx(k,1), g.ePhRx(k,1)] = getAntennaTowardAntenna(scenario.xyzTx0, scenario.xyzRx0, txSurf0, rxSurf0);
    end
    
    % Signs for summing of integrals
    iteration.signs = 1 - 2*[s>0 s<1 t>0 t<1];
    
    % Excess path lengths for Fresnel integrals
    iteration.deltaR(:,1)  = max(0,sqrt(sum(abs2(iteration.xyzEdges-scenario.xyzTx0),2)) + sqrt(sum(abs2(iteration.xyzEdges-scenario.xyzRx0),2)) - R0reflex);
        
    rMiddle = norm(iteration.xyzRP(1,:)-scenario.xyzTx0) + norm(iteration.xyzRP(1,:)-scenario.xyzRx0);
    scenario.rAmp(k,1) = norm(iteration.xyzAmp-scenario.xyzTx0) + norm(iteration.xyzAmp-scenario.xyzRx0);
    iteration.R1 = norm(iteration.xyz0(1,:)-scenario.xyzTx0);
    iteration.R2 = norm(iteration.xyz0(1,:)-scenario.xyzRx0);
    
    % Green's function times "Friis factor"
    iteration.gff    = exp(-1j*scenario.k0*rMiddle)/rMiddle * scenario.lambda/(4*pi);
    
    % Calculate the different models
    if any(strcmp(models.List, '3DFresnel'))
        gammaFresnelSingleReflection(k,1) = ourFresnelModel(scenario, iteration, 1);
        gammaFresnel(k,1) = scenario.gamma(1)*gammaFresnelSingleReflection(k,1);
    end
    if any(strcmp(models.List, 'erf'))
        gammaErfSingleReflection(k,1)     = ourErfModel(scenario, iteration, 1);
        gammaErf(k,1) = scenario.gamma(1)*gammaErfSingleReflection(k,1);
    end
    if any(strcmp(models.List, 'M-METIS'))
        gammaMetisSingleReflection(k,1)   = metisModel(scenario, iteration, 1);
        gammaMetis(k,1) = scenario.gamma(1)*gammaMetisSingleReflection(k,1);
    end
    if any(strcmp(models.List, 'METISRCS'))
        gammaMetisRCSSingleReflection(k,1)= metisRCSmodel(scenario, iteration, 1);
        gammaMetisRCS(k,1) = scenario.gamma(1)*gammaMetisRCSSingleReflection(k,1);
    end
    if any(strcmp(models.List, 'M-mmMAGIC'))
        gammaMMMagicSingleReflection(k,1) = mmMAGICModel(scenario, iteration, k, 1);
        gammaMMMagic(k,1) = scenario.gamma(1)*gammaMMMagicSingleReflection(k,1);
    end
    if any(strcmp(models.List, 'M-ITUFresnel'))
        gammaITUFresSingleReflection(k,1) = ITUfresnelModel(scenario, iteration, 1);
        gammaITUFres(k,1) = scenario.gamma(1)*gammaITUFresSingleReflection(k,1);
    end
    if any(strcmp(models.List, 'M-ITUEmpirical'))
        gammaITUEmpSingleReflection(k,1)  = ITUsemiEmpModel(scenario, iteration, k, 1);
        gammaITUEmp(k,1) = scenario.gamma(1)*gammaITUEmpSingleReflection(k,1);
    end
    
    counter = 1;
    if any(strcmp(models.List, '3DFresnel'))
        modelsTmp(counter, k) = gammaFresnel(k,1);
        counter = counter + 1;
    end
    if any(strcmp(models.List, 'erf'))
        modelsTmp(counter, k) = gammaErf(k, 1);
        counter = counter + 1;
    end
    if any(strcmp(models.List, 'M-METIS'))
        modelsTmp(counter, k) = gammaMetis(k, 1);
        counter = counter + 1;
    end
    if any(strcmp(models.List, 'METISRCS'))
        modelsTmp(counter, k) = gammaMetisRCS(k, 1)/abs(iteration.gff);
        counter = counter + 1;
    end
    if any(strcmp(models.List, 'M-mmMAGIC'))
        modelsTmp(counter, k) = gammaMMMagic(k, 1);
        counter = counter + 1;
    end
    if any(strcmp(models.List, 'M-ITUFresnel'))
        modelsTmp(counter, k) = gammaITUFres(k, 1);
        counter = counter + 1;
    end
    if any(strcmp(models.List, 'M-ITUEmpirical'))
        modelsTmp(counter, k) = gammaITUEmp(k, 1);
    end
    
    modelsTh(:, k) = modelsTmp(:, k)*...
        iteration.eThTx*iteration.eThRx*iteration.gff;
    modelsPh(:, k) = modelsTmp(:, k)*...
        iteration.ePhTx*iteration.ePhRx*iteration.gff;
end

counter = 1;
if any(strcmp(models.List, '3DFresnel'))
    models.resultsTh(:, counter) = modelsTh(counter, :);
    models.resultsPh(:, counter) = modelsPh(counter, :);
    
    counter = counter + 1;
    models.lineStyles = '- ';
end
if any(strcmp(models.List, 'erf'))
    models.resultsTh(:, counter) = modelsTh(counter, :);
    models.resultsPh(:, counter) = modelsPh(counter, :);
    
    counter = counter + 1;
    models.lineStyles = [models.lineStyles; '- '];
end
if any(strcmp(models.List, 'M-METIS'))
    models.resultsTh(:, counter) = modelsTh(counter, :);
    models.resultsPh(:, counter) = modelsPh(counter, :);
    
    counter = counter + 1;
    models.lineStyles = [models.lineStyles; '- '];
end
if any(strcmp(models.List, 'METISRCS'))
    models.resultsTh(:, counter) = modelsTh(counter, :);
    models.resultsPh(:, counter) = modelsPh(counter, :);
    
    counter = counter + 1;
    models.lineStyles = [models.lineStyles; '- '];
end
if any(strcmp(models.List, 'M-mmMAGIC'))
    models.resultsTh(:, counter) = modelsTh(counter, :);
    models.resultsPh(:, counter) = modelsPh(counter, :);
    
    counter = counter + 1;
    models.lineStyles = [models.lineStyles; '- '];
end
if any(strcmp(models.List, 'M-ITUFresnel'))
    models.resultsTh(:, counter) = modelsTh(counter, :);
    models.resultsPh(:, counter) = modelsPh(counter, :);
    
    counter = counter + 1;
    models.lineStyles = [models.lineStyles; '- '];
end
if any(strcmp(models.List, 'M-ITUEmpirical'))
    models.resultsTh(:, counter) = modelsTh(counter, :);
    models.resultsPh(:, counter) = modelsPh(counter, :);
    
    models.lineStyles = [models.lineStyles; '--'];
end
if any(strcmp(models.List, 'Free space & PEC'))
    models.lineStyles = '-';
end

if scenario.number == 2 % LoS component
    scenario.rLOS = max(0,norm(scenario.xyzRx0-scenario.xyzTx0));
    gffLOS        = exp(-1j*scenario.k0*scenario.rLOS)/scenario.rLOS * scenario.lambda/(4*pi);
    
    [eThTx, ePhTx, ~, ~] = getAntennaTowardAntenna(scenario.xyzTx0, scenario.xyzRx0, txSurf0, rxSurf0);
    LOScomponentTh = eThTx*g.eThRx.*gffLOS;
    LOScomponentPh = ePhTx*g.ePhRx.*gffLOS;
    
    models.resultsTh = models.resultsTh + LOScomponentTh;
    models.resultsPh = models.resultsPh + LOScomponentPh;
end
end

function PO = getPOresults(scenario, PO)
if PO.available == 1
    POresults = loadPOdata(scenario.number).';
else
    POresults = [];
end
PO.results = POresults;
end

function MoM = getMoMresults(scenario, MoM)
if MoM.available == 1
    MoMresults = loadMoMdata(scenario.number);
else
    MoMresults = [];
end
MoM.results = MoMresults;
end

function [errorTable, extraError] = getErrorMetric(scenario, x, y, models)
% y is the MoM reference
% x is the PO reference
% xEst is the estimation

xEst = models.resultsTh;

if y.available
    errorTable = zeros(numel(models.List),1);
    [~, b] = size(xEst);   
    idx = 1:size(y.results,1);
    
    for i = 1:b
        errorTable(i,1) = 10*log10(NMSE(y.results(idx), xEst(idx, i)));
    end
    metric  = 'momNMSEdB';    
    
    format bank
    T = table(models.List, errorTable,  'VariableNames', {'Model', metric, });
    disp(T)
    format short
    
elseif x.available
    errorTable = zeros(numel(models.List),1);
    [~, b] = size(xEst);
    
    switch scenario.number
        case {999, 9991, 7}
            w = scenario.R.^2;
            w = w/norm(w);
            for i = 1:b
                errorTable(i,1) = 10*log10(WNMSE(x.results, xEst(:, i), w));
            end
            metric       = 'poWNMSEdB';
        otherwise
            idx = 1:size(x.results,1);
            for i = 1:b
                errorTable(i,1) = 10*log10(NMSE(x.results(idx), xEst(idx, i)));
            end
            metric       = 'poNMSEdB';
    end
    
    format bank 
    T = table(models.List, errorTable, 'VariableNames', {'Model', metric});
    disp(T)
    format short
else
    errorTable = [];
    extraError = [];
end
end

function NMSE = NMSE(x, xEst)
NMSE  = sum(abs(x - xEst).^2)/sum(abs(x).^2);
end

function WNMSE = WNMSE(x, xEst, w)
WNMSE = sum(w.*abs(x - xEst).^2)/sum(w.*abs(x).^2);
end


function [scenario, PO, MoM, models] = plotResults(scenario, PO, MoM, models)
if any(scenario.number == [2, 4, 5])
    scenario.sweepVals = scenario.sweepVals/pi*180;
end

% Amplitude plot
[models.tikzAmp, PO.tikzAmp, MoM.tikzAmp]       = plotAmplitudResults(scenario, models, PO, MoM);
% Phase plot
[models.tikzPhase, PO.tikzPhase, MoM.tikzPhase] = plotPhaseResults(scenario, models, PO, MoM);

if PO.available
    models.List = [models.List; 'PO'];
end
if MoM.available
    models.List = [models.List; 'MoM'];
end
figure(999);  legend(models.List, 'Location', 'northeastoutside')
figure(1000); legend(models.List, 'Location', 'northeastoutside')
end

function [modelsTikzAmp, POtikzAmp, MoMtikzAmp] = plotAmplitudResults(scenario, models, PO, MoM)
results = models.resultsTh;

xValue = scenario.sweepVals;
modelsTikzAmp = zeros(size(results));
index = 1:numel(xValue);

figure(999);
clf

for iii = 1:numel(models.List)
    if  scenario.number == 7 % Distance attenuation removed for increasing-distance scenario
        modelsTikzAmp(:, iii) = 10*log10(abs2(results(index,iii))./scenario.R.^2);
    else
        modelsTikzAmp(:, iii) = 10*log10(abs2(results(index,iii)));
    end
    plot(xValue, modelsTikzAmp(:, iii), 'LineStyle', models.lineStyles(iii,:), 'LineWidth', 1); hold on;
end

if PO.available
    if scenario.number == 7
        POtikzAmp = 10*log10(abs2(PO.results(index,1))./scenario.R(index,1).^2);
    else
        POtikzAmp = 10*log10(abs2(PO.results(index,1)));
    end
    plot(xValue, POtikzAmp, 'k--', 'LineWidth', 1); hold on;
else
    POtikzAmp = [];
end

if MoM.available
    if scenario.number == 7
        MoMtikzAmp = 10*log10(abs2(MoM.results(index,1))./scenario.R(index,1).^2);
    else
        MoMtikzAmp = 10*log10(abs2(MoM.results(index,1)));
    end
    plot(xValue, MoMtikzAmp, '--', 'LineWidth', 1); hold on;
else
    MoMtikzAmp = [];
end

if scenario.number == 7 % Set logarithmic axes for increasing-distance scenario
    set(gca,'XScale','log')
end
getProperAxisLimitsAmplitude(scenario)
getProperXaxisLabel(scenario)
ylabel('Magnitude (dB)')
grid on
end

function [modelsTikzPhase, POtikzPhase, MoMtikzPhase] = plotPhaseResults(scenario, models, PO, MoM)

results = models.resultsTh;

xValue = scenario.sweepVals;
modelsTikzPhase = zeros(size(results));

index = 1:numel(xValue);

figure(1000);
clf
for iii = 1:numel(models.List)
    switch scenario.number
        case {7}
            phase_offset = exp(-1j*scenario.k0.*scenario.R);
        case 4
            phase_offset = exp(-1j*pi/12)*ones(numel(results));
        otherwise
            phase_offset = 1*ones(length(index));
    end
    modelsTikzPhase(:, iii) = wrapTo360(angle(results(index,iii)./phase_offset(index,1))*180/pi);
    
    plot(xValue, modelsTikzPhase(:, iii), 'LineStyle', models.lineStyles(iii,:), 'LineWidth', 1); hold on;
end

if PO.available
    if scenario.number == 4
        phase_offset = exp(-1j*pi/12)*ones(numel(results));
    end
    POtikzPhase = wrapTo360(angle(PO.results(index,1)./phase_offset(index,1))*180/pi);
    plot(xValue, POtikzPhase, 'k--', 'LineWidth', 1); hold on;
else
    POtikzPhase = [];
end

if MoM.available
    MoMtikzPhase = wrapTo360(angle(MoM.results(index,1)./phase_offset(index,1))*180/pi);
    plot(xValue, MoMtikzPhase, '--', 'LineWidth', 1); hold on;
else
    MoMtikzPhase = [];
end

if scenario.number == 7 % Set logarithmic axes for increasing-distance scenario
    set(gca,'XScale','log')
end
getProperAxisLimitsPhase(scenario)
getProperXaxisLabel(scenario);
ylabel('Phase (deg)')
grid on
end

function getProperXaxisLabel(scenario)
switch scenario.number
    case 2
        xlabel('Rx antenna pointing direction (deg)')
    case {4, 5}
        xlabel('Rectangle rotation (deg)')
    case 6
        xlabel('Closet midpoint x-offset from specular reflection point (m)')    
    case 7
        xlabel('Distance Rx antenna to rectangle centroid (m)')
end
end

function getProperAxisLimitsAmplitude(scenario)
switch scenario.number
    case 2
        xlim([-90 180]);
        ylim([-110 -45]);
    case 4
        xlim([0 180]);
        ylim([-70 -35]);
    case 5
        xlim([0 180]);        
        ylim([-60 -35]);
    case 6
        xlim([-1 1]);
        ylim([-85 -58]);    
    case 7
        ylim([-200 0]);         
    otherwise
end
end

function getProperAxisLimitsPhase(scenario)
switch scenario.number
    case 2
        xlim([-90 180]);
        ylim([0 200]);
    case {4, 5}
        xlim([0 180]);
        ylim([0 360]);
    case 6
        xlim([-1 1]);
        ylim([0 360]); 
    case 7
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
% getSpecularReflectionPoint  Image source in plane containing
% 'xyzOnSurface' with normal 'n'.
xyzTxImage = xyzTx + 2*((xyzOnSurface-xyzTx)*n.')*n; % Image of Tx
end

function [xyzRP, xyzTxImage]  = getSpecularReflectionPoint(xyzTx, xyzRx, xyzOnSurface, n)
% getSpecularReflectionPoint  Specular reflection point on plane containing
% 'xyzOnSurface' with normal 'n'.

xyzTxImage = xyzTx + 2*((xyzOnSurface-xyzTx)*n.')*n; % Image of Tx
u         = (xyzRx-xyzTxImage)/norm(xyzRx-xyzTxImage);
w         = xyzTxImage-xyzOnSurface;
s         = -(n*w.')/(n*u.');
xyzRP     = xyzTxImage + s*u;
end

function [xyz0,s0,t0] = getPointOnSurfaceClosestToRP(xyzRp,xyz1moved,t1,t2)
u  = t1/norm(t1);                       % Unit vector along first edge from first to second corner
v  = t2/norm(t2);                       % Unit vector along second edge from first to third corner
s0 = (xyzRp-xyz1moved)*u.'/norm(t1);    % Projection of vector from rectangle center to reflection point on 'u'
t0 = (xyzRp-xyz1moved)*v.'/norm(t2);    % Projection of vector from rectangle center to reflection point on 'v'
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

function res = rectangleForRotation(nSamples, freq)
surfaceWidth  = 0.3;
surfaceHeight = 0.5;
surfaceOffset = [0.25 0 -0.5*surfaceHeight];
surfScatt     = createRectangle(surfaceWidth,surfaceHeight,nSamples,freq);
surfScatt     = rotateObject(surfScatt, [1 0 0], pi/2); 
surfScatt     = moveObject(surfScatt,surfaceOffset);
res           = {surfScatt};
end

function res = rectangleForMoMorDistance(nSamples, freq)
surfaceWidth  = 0.3;
surfaceHeight = 0.5;
surfaceOffset = [0 0 0];
surfScatt     = createRectangle(surfaceWidth,surfaceHeight,nSamples,freq);
surfScatt     = rotateObject(surfScatt, [1 0 0], pi/2); 
surfScatt     = moveObject(surfScatt,surfaceOffset);
res           = {surfScatt};
end


function res = rectangleLeuven(nSamples, freq)
surfaceWidth  = 1.01;
surfaceHeight = 2;
surfaceOffset = [0 0 0.5*surfaceHeight];
surfScatt     = createRectangle(surfaceWidth,surfaceHeight,nSamples,freq);
surfScatt     = rotateObject(surfScatt, [1 0 0], pi/2); 
surfScatt     = moveObject(surfScatt,surfaceOffset);
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
% 3D Fresnel model
% [1] A. Lahuerta-Lavieja, M. Johansson, U. Gustavsson, T. A. H. Bressner, and G. A. E. Vandenbosch, 
% "Computationally-efficient millimeter-wave back-scattering models," to be published in IEEE Trans. Antennas Propag., 2020.

lambda = scenario.lambda;
deltaR = iteration.deltaR(:,m);
signs  = iteration.signs;
fresnelInts = fresnelC(sqrt(4/lambda*deltaR)) - 1j*fresnelS(sqrt(4/lambda*deltaR));
gammaFresnel = 0.5*1j*(((fresnelInts(1)*signs(1)+fresnelInts(2)*signs(2))*(fresnelInts(3)*signs(3)+fresnelInts(4)*signs(4))));
end

function gammaErf = ourErfModel(scenario, iteration, m)
% erf model [1]
% [1] A. Lahuerta-Lavieja, M. Johansson, U. Gustavsson, T. A. H. Bressner, and G. A. E. Vandenbosch, 
% "Computationally-efficient millimeter-wave back-scattering models," to be published in IEEE Trans. Antennas Propag., 2020.

lambda = scenario.lambda;
deltaR = iteration.deltaR(:,m);
signs  = iteration.signs;
erfInts = erf(sqrt(4/lambda*deltaR));
gammaErf = ((erfInts(1)*signs(1)+erfInts(2)*signs(2))*(erfInts(3)*signs(3)+erfInts(4)*signs(4)))*0.25;
end

function gammaMetis = metisModel(scenario, iteration, m)
% M-METIS model
% Based on METIS D1.4 Section C.1.4: Shadowing objects
% url: https://www.metis2020.com/wp-content/uploads/METIS_D1.4_v3.pdf

lambda = scenario.lambda;
deltaR = iteration.deltaR(:,m);
signs  = iteration.signs;
metisTerms = atan(0.5*pi*sqrt(pi/lambda*deltaR))/pi;
gammaMetis = ((metisTerms(1)*signs(1)+metisTerms(2)*signs(2))*(metisTerms(3)*signs(3)+metisTerms(4)*signs(4)));
end

function [gammaMetisRCS] =  metisRCSmodel(scenario, iteration, m)
% METIS RCS scattering model
% Based on METIS D1.4 Section C.1.5 Scattering objects
% url: https://www.metis2020.com/wp-content/uploads/METIS_D1.4_v3.pdf

lambda = scenario.lambda;
deltaR = iteration.deltaR(:,m);
signs  = iteration.signs;
R1     = iteration.R1;
R2     = iteration.R2;
area   = iteration.area(:,m);

radius = sqrt(area/pi);
RCS    = pi*radius^2;
metisTerms = atan(0.5*pi*sqrt(pi/lambda*deltaR))/pi;
gammaMetisRCS0 = ((metisTerms(1)*signs(1)+metisTerms(2)*signs(2))*(metisTerms(3)*signs(3)+metisTerms(4)*signs(4)));
Ssc             = 1/(4*pi*(R1)^2);
correction = (4*pi*area)/lambda^2;
gammaMetisRCS =  sqrt(correction*Ssc*RCS*(lambda/(4*pi*R2))^2*(1-gammaMetisRCS0)^2);
end

function gammaMMMagic = mmMAGICModel(scenario, iteration, k, m)
% M-mmMAGIC model
% Based on mmMAGIC D2.2 Section 4.6.3 Blockage
% url: https://bscw.5g-mmmagic.eu/pub/bscw.cgi/d202656/mmMAGIC_D2-2.pdf

lambda     = scenario.lambda;
k0         = scenario.k0;
xyzRx0     = scenario.xyzRx0;
R          = scenario.R(k);
deltaR     = iteration.deltaR(:,m);
signs      = iteration.signs;
xyzEdges   = iteration.xyzEdges;
xyzTxImage = iteration.xyzTxImage(m,:);

metisTerms = atan(0.5*pi*sqrt(pi/lambda*deltaR))/pi;

edgesVectorsTx = xyzEdges - xyzTxImage;
edgesVectorsRx = xyzRx0 - xyzEdges;
cos_mmMAGIC_phi = dot(edgesVectorsTx, edgesVectorsRx, 2)./(sqrt(sum(abs2(edgesVectorsTx),2)).*sqrt(sum(abs2(edgesVectorsRx),2)));

mmMAGIC_phase_projection = exp(-1j*k0*(sqrt(sum(abs2(xyzEdges-xyzTxImage),2)) + sqrt(sum(abs2(xyzEdges-xyzRx0),2))));
mmMAGIC_phase_R = exp(-1j*k0*R);
mmMAGIC_F = cos_mmMAGIC_phi.*(0.5-metisTerms);
mmMAGICTerms = (0.5 - mmMAGIC_phase_projection/mmMAGIC_phase_R.*mmMAGIC_F);

gammaMMMagic = ((mmMAGICTerms(1)*signs(1)+mmMAGICTerms(2)*signs(2))*(mmMAGICTerms(3)*signs(3)+mmMAGICTerms(4)*signs(4)));
end

function gammaITUFres = ITUfresnelModel(scenario, iteration, m)
% Based on ITU Recommendation P.526-14 Section 5.2.1.1: Fresnel integral method
% url: https://www.itu.int/dms_pubrec/itu-r/rec/p/R-REC-P.526-14-201801-I!!PDF-E.pdf

lambda     = scenario.lambda;
xyzRx0     = scenario.xyzRx0;
xyzEdges   = iteration.xyzEdges;
xyzTxImage = iteration.xyzTxImage(m,:);
xyzRP      = iteration.xyzRP(m,:);
xyzTx0     = xyzTxImage;

t1 = xyzEdges(2,:) - xyzEdges(1,:);
t1 = t1/norm(t1);
t2 = xyzEdges(4,:) - xyzEdges(3,:);
t2 = t2/norm(t2);
n  = cross(t1,t2);
n  = n/norm(n);

zt = dot(xyzRP-xyzTx0, n);
zr = dot(xyzRx0-xyzRP, n);
xy = xyzEdges - xyzRP;

rho_r1 = sum(xyzEdges(1:2,:) .* t1, 2) - sum(xyzRx0.* t1, 2);
rho_t1 = sum(xyzEdges(1:2,:) .* t1, 2) - sum(xyzTx0.* t1, 2);
rho_r2 = sum(xyzEdges(3:4,:) .* t2, 2) - sum(xyzRx0.* t2, 2);
rho_t2 = sum(xyzEdges(3:4,:) .* t2, 2) - sum(xyzTx0.* t2, 2);

phi = zeros(4,1);
phi(1:2) = atan2(rho_r1, zr) + atan2(rho_t1, zt);
phi(3:4) = atan2(rho_r2, zr) + atan2(rho_t2, zt);

x = dot(xy, t1(ones(size(xy,1),1),:), 2);
y = dot(xy, t2(ones(size(xy,1),1),:), 2);
xy = [x(1:2); y(3:4)];

v_ITU = sign(xy).*sqrt(2/lambda.*abs(xy).^1.18*abs(1/zr + 1/zt)^0.18.*abs(phi).^0.82);

Cx = fresnelC(v_ITU(2)) - fresnelC(v_ITU(1));
Cy = fresnelC(v_ITU(4)) - fresnelC(v_ITU(3));
Sx = fresnelS(v_ITU(2)) - fresnelS(v_ITU(1));
Sy = fresnelS(v_ITU(4)) - fresnelS(v_ITU(3));

gammaITUFres = conj((Cx*Sy + Sx*Cy) + 1j*(Sx*Sy - Cx*Cy))*0.5;
end

function gammaITUEmp = ITUsemiEmpModel(scenario, iteration, k, m)
% Based on ITU Recommendation P.526-14 Section 5.2.1.2: Semi-empirical method
% url: https://www.itu.int/dms_pubrec/itu-r/rec/p/R-REC-P.526-14-201801-I!!PDF-E.pdf

lambda     = scenario.lambda;
xyzRx0     = scenario.xyzRx0;
R          = scenario.R(k);
deltaR     = iteration.deltaR(:,m);
xyzEdges   = iteration.xyzEdges;
xyzTxImage = iteration.xyzTxImage(m,:);
xyzRP      = iteration.xyzRP(m,:);

xyzTx0 = xyzTxImage;

t1 = xyzEdges(2,:) - xyzEdges(1,:);
t1 = t1/norm(t1);
t2 = xyzEdges(4,:) - xyzEdges(3,:);
t2 = t2/norm(t2);
n  = cross(t1,t2);
n  = n/norm(n);

zt = dot(xyzRP-xyzTx0, n);
zr = dot(xyzRx0-xyzRP, n);
xy = xyzEdges - xyzRP;

rho_r1 = sum(xyzEdges(1:2,:) .* t1, 2) - sum(xyzRx0.* t1, 2);
rho_t1 = sum(xyzEdges(1:2,:) .* t1, 2) - sum(xyzTx0.* t1, 2);
rho_r2 = sum(xyzEdges(3:4,:) .* t2, 2) - sum(xyzRx0.* t2, 2);
rho_t2 = sum(xyzEdges(3:4,:) .* t2, 2) - sum(xyzTx0.* t2, 2);

phi = zeros(4,1);
phi(1:2) = atan2(rho_r1, zr) + atan2(rho_t1, zt);
phi(3:4) = atan2(rho_r2, zr) + atan2(rho_t2, zt);

x = dot(xy, t1(ones(size(xy,1),1),:), 2);
y = dot(xy, t2(ones(size(xy,1),1),:), 2);
xy = [x(1:2); y(3:4)];

k0 = 2*pi/lambda;
metisTerms = atan(2*1.4*sqrt(deltaR/lambda));

ITU_phase_projection = exp(-1j*k0*(sqrt(sum(abs2(xyzEdges-xyzTxImage),2)) + sqrt(sum(abs2(xyzEdges-xyzRx0),2))));
ITU_phase_R = exp(-1j*k0*R);
ITU_G = cos(phi*0.5).*(0.5-metisTerms/pi);
ITUTerms = (0.5 - ITU_phase_projection/ITU_phase_R.*ITU_G);

gammaITUEmp = ((ITUTerms(1)*sign(xy(1))-ITUTerms(2)*sign(xy(2)))*(ITUTerms(3)*sign(xy(3))-ITUTerms(4)*sign(xy(4))));
end

function result = abs2(x)
result = real(x).^2 + imag(x).^2;
end

function [surface, nVerts] = createRectangle(width,height,nVerts, frequency)
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

function obj = moveObject(obj, dist) 
switch obj.type 
    case 'antenna'        
        obj.geomT = {dist obj.geomT{:}};
    otherwise 
        obj.vertices = bsxfun(@plus, obj.vertices, dist);         
        obj.centroid = bsxfun(@plus, obj.centroid, dist);
end 
end

function obj = rotateObject(obj, rotAxis, alpha) 
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
% NardaV637_TAP
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
