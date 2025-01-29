%%%%% PVK vs Si panel temperature study     %%%%%
% Solar Energy Materials and Systems, University of Turku

% Assume 0.9 total absorption for energies above band gap (0.1 total reflectance)
% Assume 0.2 total absorption for energies below band gap (0.8 total reflectance)
% Assume 20% efficiency for both cells

% Add comsol modelling

errorCountLimit=5;
% PEROVSKITE

% Parameters
% Natural constants
consts = ([]);
consts.('h') = 6.626*10^-34;    % Js
consts.('h_eV') = 4.136*10^-15; % eVs
consts.('c') = 2.998*10^8;      % m/s
consts.('kB') = 1.381*10^-23;   % J/K
consts.('q') = 1.602*10^-19;    % C
consts.('SB') = 5.670*10^-8;    % W/(m^2K^4);

% Band gap
% MAPbI3
EgRef = 1.57;
EgRefT = 298.15; % K
EgShiftVsT = 0.25*1e-3; % eV/K, https://doi.org/10.1021/acs.jpclett.6b01207

% Total efficiencies
EffRefT = 298.15; % K
EffRefG = 1000; % W/m^2
VocRefG = 1.1; % W
n=1.3;

% % PSC
effRef = 0.2;
EffShiftVsTPVK = -0.17; % %/K

% Spectral absorption coeff above and below the respective band gaps
absCoeffAboveGap = 0.95;
absCoeffBelowGap = 0.2;

% Irradiation data
AM15 = load('AM15.mat').AM15;

% Modelled wavelengths
dLambda = 1; % nm
lambda = 300:dLambda:4000; 
lambda = transpose(lambda);

% Modelled spectrum
lambda2 = AM15{:,1};                        % Wavelength (nm)
B = AM15{:,2};                              % Radiance (Wm^-2nm^-1)
Bsolar = spline(lambda2,B,lambda);
BsolarTot = sum(Bsolar*dLambda);

% Folder
pathFolder ='.\PSCAllCond';

% Error log
FolderErrorLog = [pathFolder,'\errorLog'];
mkdir(FolderErrorLog);

% Initial temperature 
Tstc = 293.15; % K
Tcell = Tstc; % K

% Weather conditions (simultaneously changing)
Tamb = [0 10 20 30 40]; % C
Irradiance = [200 400 600 800 1000 1200]; % W/m^2
v_wind = [0 2 4 6 8 10]; % m/s (unitless in comsol model)

% Panel design
cellDistance = 1; % cm
Acells = 72*0.156^2; % m^2

% Temperature convergence
TendCriterion = 0.1; % K
countTlimit = 10;

% Result storage
% % Cell and module temperatures 

TaveCells_vs_PSCAll = zeros(length(Tamb),length(Irradiance),length(v_wind),length(cellDistance),countTlimit); % Average temperature of PSC cells
TcellAveTop_vs_PSCAll = zeros(length(Tamb),length(Irradiance),length(v_wind),length(cellDistance),countTlimit); % Average temperature of the top surface of PSC module
TcellAveBack_vs_PSCAll = zeros(length(Tamb),length(Irradiance),length(v_wind),length(cellDistance),countTlimit); % Average temperature of the back surface of PSC module
TaveMod_vs_PSCAll = zeros(length(Tamb),length(Irradiance),length(v_wind),length(cellDistance),countTlimit); % Average temperature of PSC module (all domains)
TmaxCells_vs_PSCAll = zeros(length(Tamb),length(Irradiance),length(v_wind),length(cellDistance),countTlimit); % Maximum temperature of PSC cells

% % Efficiency
Efficiency_vs_PSCAll = zeros(length(Tamb),length(Irradiance),length(v_wind),length(cellDistance),countTlimit); 

% % Cutlines
CutLineCellsResults_vs_PSCAll = cell(length(Tamb),length(Irradiance),length(v_wind),length(cellDistance),countTlimit); % Temperature through the module at the level of the cells
CutLineTopResults_vs_PSCAll = cell(length(Tamb),length(Irradiance),length(v_wind),length(cellDistance),countTlimit); % Temperature through the module at the level of the module top surface
CutLineBackResults_vs_PSCAll = cell(length(Tamb),length(Irradiance),length(v_wind),length(cellDistance),countTlimit); % Temperature through the module at the level of the module back surface
CutLineCellsXResults_vs_PSCAll = cell(length(Tamb),length(Irradiance),length(v_wind),length(cellDistance),countTlimit); % The points at which the "CutLineCellsResults_vs_PSCAll" temperatures are determined
CutLineTopXResults_vs_PSCAll = cell(length(Tamb),length(Irradiance),length(v_wind),length(cellDistance),countTlimit); % The points at which the "CutLineTopResults_vs_PSCAll" temperatures are determined
CutLineBackXResults_vs_PSCAll = cell(length(Tamb),length(Irradiance),length(v_wind),length(cellDistance),countTlimit); % The points at which the "CutLineBackResults_vs_PSCAll" temperatures are determined

% % Absorption, Heat Production, Convective and radiative heat flux
AbsorptionCoeff_vs_PSCAll = zeros(length(Tamb),length(Irradiance),length(v_wind),length(cellDistance),countTlimit);
HeatProduction_vs_PSCAll = zeros(length(Tamb),length(Irradiance),length(v_wind),length(cellDistance),countTlimit);
ConvAveTop_vs_PSCAll = zeros(length(Tamb),length(Irradiance),length(v_wind),length(cellDistance),countTlimit);
ConvAveBack_vs_PSCAll = zeros(length(Tamb),length(Irradiance),length(v_wind),length(cellDistance),countTlimit);
ConvAveSide_vs_PSCAll = zeros(length(Tamb),length(Irradiance),length(v_wind),length(cellDistance),countTlimit);
RadAveTop_vs_PSCAll = zeros(length(Tamb),length(Irradiance),length(v_wind),length(cellDistance),countTlimit);
RadAveBack_vs_PSCAll = zeros(length(Tamb),length(Irradiance),length(v_wind),length(cellDistance),countTlimit);
RadAveSide_vs_PSCAll = zeros(length(Tamb),length(Irradiance),length(v_wind),length(cellDistance),countTlimit);

QcellTot_vs_PSCAll = zeros(length(Tamb),length(Irradiance),length(v_wind),length(cellDistance),countTlimit);

% Loop over environmental conditions
for idxTamb = 1:length(Tamb)
for idxIrr = 1:length(Irradiance)
for idxv_wind = 1:length(v_wind)
for idxcellDistance = 1:length(cellDistance)

%%% Loop over temperature iteration
countT = 0;
Tchange = 999;
errorCount = 0;
Tcell = Tstc; % K

while (Tchange > TendCriterion && countT < countTlimit)
% Start time
tStartRound = tic;
try 
    % Band gap at a temperature
    EgAtT = EgRef+EgShiftVsT*(Tcell-EgRefT);
    EgWl = consts.h_eV*consts.c*10^9/EgAtT; % nm

    % Efficiency at current temperature

    % Absorption
    A = ones(size(Bsolar))*absCoeffBelowGap;
    A(lambda<EgWl) = absCoeffAboveGap;

    % Total absorption coefficient
    % Spectrum weighted average of normalized wavelength specific absorption
    totAbsCoeff = sum(Bsolar.*A*dLambda)/BsolarTot;

    % Electrical efficiency
    EffRefAtT = effRef*(1+(n*consts.kB*EffRefT)/(VocRefG*consts.q)*log(Irradiance(idxIrr)/EffRefG))*(1+EffShiftVsTPVK*(Tcell-EffRefT)/100);
% 
    % EffRefAtT = effRef*(1+EffShiftVsTPVK*(Tcell-EffRefT)/100);
    effEref = EffRefAtT/totAbsCoeff;

    % Heat production
    QprodCoeff = (1-effEref)*totAbsCoeff;

    % Irradiance effect
    Qcell = Irradiance(idxIrr)*QprodCoeff; % W/m^2

    % Build Comsol model
    [model,TcellAve,TcellAveTop,TcellAveBack,TmodAve,TcellMax,Conv_AveTop,Conv_AveBack,Conv_AveSide,Rad_AveTop,Rad_AveBack,Rad_AveSide,CutLineCells,CutLineTop,CutLineBack,CutLineCellsX,CutLineTopX,CutLineBackX] = PSCpanelTcomsolModel_v8(Tamb(idxTamb),v_wind(idxv_wind),cellDistance(idxcellDistance),Qcell,pathFolder,Tcell);

    % Extract results
    TcellUpdate = TcellAve;

    % Temperature convergence?
    Tchange = abs(TcellUpdate-Tcell);

    % Update cell temperature
    Tcell = TcellUpdate;  % Cell average

    % T iteration limit reached
    countT = countT+1;
    if countT >= countTlimit
        warning('Tcell iteration limit reached.') % TO DO: save into log file
    end

    % Save results
    TaveCells_vs_PSCAll(idxTamb,idxIrr,idxv_wind,idxcellDistance,countT) = Tcell;
    TcellAveTop_vs_PSCAll(idxTamb,idxIrr,idxv_wind,idxcellDistance,countT) = TcellAveTop;
    TcellAveBack_vs_PSCAll(idxTamb,idxIrr,idxv_wind,idxcellDistance,countT) = TcellAveBack;
    TaveMod_vs_PSCAll(idxTamb,idxIrr,idxv_wind,idxcellDistance,countT) = TmodAve;
    TmaxCells_vs_PSCAll(idxTamb,idxIrr,idxv_wind,idxcellDistance,countT) = TcellMax;
    %
    Efficiency_vs_PSCAll(idxTamb,idxIrr,idxv_wind,idxcellDistance,countT) = EffRefAtT;
    %
    CutLineCellsResults_vs_PSCAll{idxTamb,idxIrr,idxv_wind,idxcellDistance,countT} = CutLineCells;
    CutLineTopResults_vs_PSCAll{idxTamb,idxIrr,idxv_wind,idxcellDistance,countT} = CutLineTop;
    CutLineBackResults_vs_PSCAll{idxTamb,idxIrr,idxv_wind,idxcellDistance,countT} = CutLineBack;
    CutLineCellsXResults_vs_PSCAll{idxTamb,idxIrr,idxv_wind,idxcellDistance,countT} = CutLineCellsX;
    CutLineTopXResults_vs_PSCAll{idxTamb,idxIrr,idxv_wind,idxcellDistance,countT} = CutLineTopX;
    CutLineBackXResults_vs_PSCAll{idxTamb,idxIrr,idxv_wind,idxcellDistance,countT} = CutLineBackX;
    %
    AbsorptionCoeff_vs_PSCAll(idxTamb,idxIrr,idxv_wind,idxcellDistance,countT) = totAbsCoeff;
    HeatProduction_vs_PSCAll(idxTamb,idxIrr,idxv_wind,idxcellDistance,countT) = Qcell;
    ConvAveTop_vs_PSCAll(idxTamb,idxIrr,idxv_wind,idxcellDistance,countT) = Conv_AveTop;
    ConvAveBack_vs_PSCAll(idxTamb,idxIrr,idxv_wind,idxcellDistance,countT) = Conv_AveBack;
    ConvAveSide_vs_PSCAll(idxTamb,idxIrr,idxv_wind,idxcellDistance,countT) = Conv_AveSide;
    RadAveTop_vs_PSCAll(idxTamb,idxIrr,idxv_wind,idxcellDistance,countT) = Rad_AveTop;
    RadAveBack_vs_PSCAll(idxTamb,idxIrr,idxv_wind,idxcellDistance,countT) = Rad_AveBack;
    RadAveSide_vs_PSCAll(idxTamb,idxIrr,idxv_wind,idxcellDistance,countT) = Rad_AveSide;

    % Total heat production
    QcellTot_vs_PSCAll(idxTamb,idxIrr,idxv_wind,idxcellDistance,countT) = Qcell*Acells; % W
    errorCount = 0;

catch ME
    errorCount = errorCount+1;
    % Write error to a file
    fid = fopen([FolderErrorLog,'\errorCount_',num2str(errorCount),'.txt'],'w');
    fprintf(fid,'%s',ME.getReport('extended','hyperlinks','off'));
    fclose(fid);

    pause(5)
    if errorCount>errorCountLimit
        keyboard
    end
end
% Progress report
tElapsedRound = toc(tStartRound);
disp(['Runtime, last round: ',num2str(tElapsedRound/60),' min.'])
end
end
end
end
end

varsRemove = {'ME','model'};
clear(varsRemove{:});

% Heat conservation
Qout = ConvAveTop_vs_PSCAll+ConvAveBack_vs_PSCAll+ConvAveSide_vs_PSCAll+RadAveTop_vs_PSCAll+RadAveBack_vs_PSCAll+RadAveSide_vs_PSCAll;
Qdiff = QcellTot_vs_PSCAll+Qout; % plus due to the signs
disp(['Maximum difference between Qin and Qout: ',num2str(max(Qdiff./QcellTot_vs_PSCAll*100)),'%.'])

% Save workspace
save([pathFolder,'\workspace.mat'])


% SILICON

%%% Parameters
% Natural constants
consts = ([]);
consts.('h') = 6.626*10^-34;    % Js
consts.('h_eV') = 4.136*10^-15; % eVs
consts.('c') = 2.998*10^8;      % m/s
consts.('kB') = 1.381*10^-23;   % J/K
consts.('q') = 1.602*10^-19;    % C
consts.('SB') = 5.670*10^-8;    % W/(m^2K^4);

% Band gap
% Si
EgRef = 1.206; % eV, https://doi.org/10.1063/1.345414
EgRefT = 0; % K, https://doi.org/10.1063/1.345414
EgShiftVsT = -0.273*1e-3; % eV/K, https://doi.org/10.1063/1.345414

% Si
effRef = 0.149; % 0.2;
EffShiftVsTSi = -0.41; % -0.35; % %/K

% % Total efficiencies
EffRefT = 298.15; % K
EffRefG = 1000; % W/m^2
VocRefG = 0.6; % W
n=1.3;

% Spectral absorption coeff above and below the respective band gaps
absCoeffAboveGap = 0.95;
absCoeffBelowGap = 0.2;

% Irradiation data
AM15 = load('AM15.mat').AM15;

% Modelled wavelengths
dLambda = 1; % nm
lambda = 300:dLambda:4000; 
lambda = transpose(lambda);

% Modelled spectrum
lambda2 = AM15{:,1};                        % Wavelength (nm)
B = AM15{:,2};                              % Radiance (Wm^-2nm^-1)
Bsolar = spline(lambda2,B,lambda);
BsolarTot = sum(Bsolar*dLambda);

% Folder
pathFolder ='.\SiAllCond';

% Error log
FolderErrorLog = [pathFolder,'\errorLog'];
mkdir(FolderErrorLog);

% Initial temperature 
Tstc = 293.15; % K

% Weather conditions (simultaneously changing)
Tamb = [0 10 20 30 40]; % C
Irradiance = [200 400 600 800 1000 1200]; % W/m^2
v_wind = [0 2 4 6 8 10]; % m/s (unitless in comsol model)

cellDistance = 1; % cm

Acells = 72*0.156^2; % m^2

% Temperature convergence
TendCriterion = 0.1; % K
countTlimit = 10;

% Result storage
% % Cell and module temperatures 
TaveCells_vs_SiAll = zeros(length(Tamb),length(Irradiance),length(v_wind),length(cellDistance),countTlimit); % Average temperature of Si cells
TcellAveTop_vs_SiAll = zeros(length(Tamb),length(Irradiance),length(v_wind),length(cellDistance),countTlimit); % Average temperature of the top surface of Si module
TcellAveBack_vs_SiAll = zeros(length(Tamb),length(Irradiance),length(v_wind),length(cellDistance),countTlimit); % Average temperature of the back surface of Si module
TaveMod_vs_SiAll = zeros(length(Tamb),length(Irradiance),length(v_wind),length(cellDistance),countTlimit); % Average temperature of Si module (all domains)
TmaxCells_vs_SiAll = zeros(length(Tamb),length(Irradiance),length(v_wind),length(cellDistance),countTlimit); % Maximum temperature of Si cells

% % Efficiency
Efficiency_vs_SiAll = zeros(length(Tamb),length(Irradiance),length(v_wind),length(cellDistance),countTlimit);

% % Cutlines
CutLineCellsResults_vs_SiAll = cell(length(Tamb),length(Irradiance),length(v_wind),length(cellDistance),countTlimit); % Temperature through the module at the level of the cells
CutLineTopResults_vs_SiAll = cell(length(Tamb),length(Irradiance),length(v_wind),length(cellDistance),countTlimit); % Temperature through the module at the level of the module top surface
CutLineBackResults_vs_SiAll = cell(length(Tamb),length(Irradiance),length(v_wind),length(cellDistance),countTlimit); % Temperature through the module at the level of the module back surface
CutLineCellsXResults_vs_SiAll = cell(length(Tamb),length(Irradiance),length(v_wind),length(cellDistance),countTlimit); % The points at which the "CutLineCellsResults_vs_SiAll" temperatures are determined
CutLineTopXResults_vs_SiAll = cell(length(Tamb),length(Irradiance),length(v_wind),length(cellDistance),countTlimit); % The points at which the "CutLineTopResults_vs_SiAll" temperatures are determined
CutLineBackXResults_vs_SiAll = cell(length(Tamb),length(Irradiance),length(v_wind),length(cellDistance),countTlimit); % The points at which the "CutLineBackResults_vs_SiAll" temperatures are determined

% % Absorption, Heat Production, Convective and radiative heat flux
AbsorptionCoeff_vs_SiAll = zeros(length(Tamb),length(Irradiance),length(v_wind),length(cellDistance),countTlimit);
HeatProduction_vs_SiAll = zeros(length(Tamb),length(Irradiance),length(v_wind),length(cellDistance),countTlimit);
ConvAveTop_vs_SiAll = zeros(length(Tamb),length(Irradiance),length(v_wind),length(cellDistance),countTlimit);
ConvAveBack_vs_SiAll = zeros(length(Tamb),length(Irradiance),length(v_wind),length(cellDistance),countTlimit);
ConvAveSide_vs_SiAll = zeros(length(Tamb),length(Irradiance),length(v_wind),length(cellDistance),countTlimit);
RadAveTop_vs_SiAll = zeros(length(Tamb),length(Irradiance),length(v_wind),length(cellDistance),countTlimit);
RadAveBack_vs_SiAll = zeros(length(Tamb),length(Irradiance),length(v_wind),length(cellDistance),countTlimit);
RadAveSide_vs_SiAll = zeros(length(Tamb),length(Irradiance),length(v_wind),length(cellDistance),countTlimit);

QcellTot_vs_SiAll = zeros(length(Tamb),length(Irradiance),length(v_wind),length(cellDistance),countTlimit);

% Loop over environmental conditions
for idxTamb = 1:length(Tamb)
for idxIrr = 1:length(Irradiance)
for idxv_wind = 1:length(v_wind)
for idxcellDistance = 1:length(cellDistance)
% for idxkcoeff = 1:length(k_coeff)

%%% Loop over temperature iteration
countT = 0;
Tchange = 999;
errorCount = 0;
Tcell = Tstc; % K
while (Tchange > TendCriterion && countT < countTlimit)
% Start time
tStartRound = tic;
try 
    % Band gap at a temperature
    EgAtT = EgRef+EgShiftVsT*(Tcell-EgRefT);
    EgWl = consts.h_eV*consts.c*10^9/EgAtT; % nm

    % Absorption
    A = ones(size(Bsolar))*absCoeffBelowGap;
    A(lambda<EgWl) = absCoeffAboveGap;

    % Total absorption coefficient
    % Spectrum weighted average of normalized wavelength specific absorption
    totAbsCoeff = sum(Bsolar.*A*dLambda)/BsolarTot;

    % Electrical efficiency
    EffRefAtT = effRef*(1+(n*consts.kB*EffRefT)/(VocRefG*consts.q)*log((Irradiance(idxIrr))/(EffRefG)))*(1+EffShiftVsTSi*(Tcell-EffRefT)/100);

    % EffRefAtT = effRef*(1+EffShiftVsTSi*(Tcell-EffRefT)/100);
    effEref = EffRefAtT/totAbsCoeff;
    % effEref = effRef/totAbsCoeff;

    % Heat production
    QprodCoeff = (1-effEref)*totAbsCoeff;


    % Irradiance effect
    Qcell = Irradiance(idxIrr)*QprodCoeff; % W/m^2

    % Build Comsol model
    [model,TcellAve,TcellAveTop,TcellAveBack,TmodAve,TcellMax,Conv_AveTop,Conv_AveBack,Conv_AveSide,Rad_AveTop,Rad_AveBack,Rad_AveSide,CutLineCells,CutLineTop,CutLineBack,CutLineCellsX,CutLineTopX,CutLineBackX] = SipanelTcomsolModel_v6(Tamb(idxTamb),v_wind(idxv_wind),cellDistance(idxcellDistance),Qcell,pathFolder,Tcell);
    
    % Extract results
    TcellUpdate = TcellAve;

    % Temperature convergence?
    Tchange = abs(TcellUpdate-Tcell);

    % Update cell temperature
    Tcell = TcellUpdate;  % Cell average

    % T iteration limit reached
    countT = countT+1;
    if countT >= countTlimit
        warning('Tcell iteration limit reached.') % TO DO: save into log file
    end

    % Save results
    TaveCells_vs_SiAll(idxTamb,idxIrr,idxv_wind,idxcellDistance,countT) = Tcell;
    TcellAveTop_vs_SiAll(idxTamb,idxIrr,idxv_wind,idxcellDistance,countT) = TcellAveTop;
    TcellAveBack_vs_SiAll(idxTamb,idxIrr,idxv_wind,idxcellDistance,countT) = TcellAveBack;
    TaveMod_vs_SiAll(idxTamb,idxIrr,idxv_wind,idxcellDistance,countT) = TmodAve;
    TmaxCells_vs_SiAll(idxTamb,idxIrr,idxv_wind,idxcellDistance,countT) = TcellMax;
    %
    Efficiency_vs_SiAll(idxTamb,idxIrr,idxv_wind,idxcellDistance,countT)= EffRefAtT;
    %
    CutLineCellsResults_vs_SiAll{idxTamb,idxIrr,idxv_wind,idxcellDistance,countT} = CutLineCells;
    CutLineTopResults_vs_SiAll{idxTamb,idxIrr,idxv_wind,idxcellDistance,countT} = CutLineTop;
    CutLineBackResults_vs_SiAll{idxTamb,idxIrr,idxv_wind,idxcellDistance,countT} = CutLineBack;
    CutLineCellsXResults_vs_SiAll{idxTamb,idxIrr,idxv_wind,idxcellDistance,countT} = CutLineCellsX;
    CutLineTopXResults_vs_SiAll{idxTamb,idxIrr,idxv_wind,idxcellDistance,countT} = CutLineTopX;
    CutLineBackXResults_vs_SiAll{idxTamb,idxIrr,idxv_wind,idxcellDistance,countT} = CutLineBackX;
    %
    AbsorptionCoeff_vs_SiAll(idxTamb,idxIrr,idxv_wind,idxcellDistance,countT) = totAbsCoeff;
    HeatProduction_vs_SiAll(idxTamb,idxIrr,idxv_wind,idxcellDistance,countT) = Qcell;
    ConvAveTop_vs_SiAll(idxTamb,idxIrr,idxv_wind,idxcellDistance,countT) = Conv_AveTop;
    ConvAveBack_vs_SiAll(idxTamb,idxIrr,idxv_wind,idxcellDistance,countT) = Conv_AveBack;
    ConvAveSide_vs_SiAll(idxTamb,idxIrr,idxv_wind,idxcellDistance,countT) = Conv_AveSide;
    RadAveTop_vs_SiAll(idxTamb,idxIrr,idxv_wind,idxcellDistance,countT) = Rad_AveTop;
    RadAveBack_vs_SiAll(idxTamb,idxIrr,idxv_wind,idxcellDistance,countT) = Rad_AveBack;
    RadAveSide_vs_SiAll(idxTamb,idxIrr,idxv_wind,idxcellDistance,countT) = Rad_AveSide;

    % Total heat production
    QcellTot_vs_SiAll(idxTamb,idxIrr,idxv_wind,idxcellDistance,countT) = Qcell*Acells; % W
    errorCount = 0;

catch ME
    errorCount = errorCount+1;
    % Write error to a file
    fid = fopen([FolderErrorLog,'\errorCount_',num2str(errorCount),'.txt'],'w');
    fprintf(fid,'%s',ME.getReport('extended','hyperlinks','off'));
    fclose(fid);

    pause(5)
    if errorCount>errorCountLimit
        keyboard
    end
end
% Progress report
tElapsedRound = toc(tStartRound);
disp(['Runtime, last round: ',num2str(tElapsedRound/60),' min.'])
end
end
end
end
end
% end


varsRemove = {'ME','model'};
clear(varsRemove{:});

% Heat conservation
Qout = ConvAveTop_vs_SiAll+ConvAveBack_vs_SiAll+ConvAveSide_vs_SiAll+RadAveTop_vs_SiAll+RadAveBack_vs_SiAll+RadAveSide_vs_SiAll;
Qdiff = QcellTot_vs_SiAll+Qout; % plus due to the signs
disp(['Maximum difference between Qin and Qout: ',num2str(max(Qdiff./QcellTot_vs_SiAll*100)),'%.'])

% Save workspace
save([pathFolder,'\workspace.mat'])
