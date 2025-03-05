% PVK vs Si panel temperature study
% Solar Energy Materials and Systems, University of Turku

% Assume 0.9 total absorption for energies above band gap (0.1 total reflectance)
% Assume 0.2 total absorption for energies below band gap (0.8 total reflectance)
% Assume 20% efficiency for both cells
% Add comsol modelling

errorCountLimit=5;
% % PEROVSKITE % 

% %%% Parameters
% % Natural constants
consts = ([]);
consts.('h') = 6.626*10^-34;    % Js
consts.('h_eV') = 4.136*10^-15; % eVs
consts.('c') = 2.998*10^8;      % m/s
consts.('kB') = 1.381*10^-23;   % J/K
consts.('q') = 1.602*10^-19;    % C
consts.('SB') = 5.670*10^-8;    % W/(m^2K^4);

% % Band gap
% % MAPbI3
EgRef = 1.57;
EgRefT = 298.15; % K
EgShiftVsT = 0.25*1e-3; % eV/K, https://doi.org/10.1021/acs.jpclett.6b01207

% Total efficiencies
EffRefT = 298.15; % K
EffRefG = 1000; % W/m^2
VocRefG = 1.1; % W
n=1.3;

% PSC
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
% 

% Ambient temperature % 

% Folder
pathFolder ='.\PSCTamb';
% Error log
FolderErrorLog = [pathFolder,'\errorLog'];
mkdir(FolderErrorLog);

% Initial temperature 
Tstc = 293.15; % K
Tcell = Tstc; % K

% Weather conditions
Tamb = [-20 -10 0 10 20 30 40 50]; % C
Irradiance = 800; % W/m^2
v_wind = 1; % m/s (unitless in comsol model)

% Panel design
cellDistance = 1; %cm

% Temperature convergence
TendCriterion = 0.1; % K
countTlimit = 10;

% Cell and module temperatures 
TaveCells_vs_PSCTamb = zeros(length(Tamb),countTlimit);
TAveTop_vs_PSCTamb = zeros(length(Tamb),countTlimit);
TAveBack_vs_PSCTamb = zeros(length(Tamb),countTlimit);
TaveMod_vs_PSCTamb = zeros(length(Tamb),countTlimit);
TmaxCells_vs_PSCTamb = zeros(length(Tamb),countTlimit);

% Efficiency
Efficiency_vs_PSCTamb = zeros(length(Tamb),countTlimit);

% Cutlines
CutLineCellsResults_vs_PSCTamb = cell(length(Tamb),countTlimit);
CutLineTopResults_vs_PSCTamb = cell(length(Tamb),countTlimit);
CutLineBackResults_vs_PSCTamb = cell(length(Tamb),countTlimit);
CutLineCellsXResults_vs_PSCTamb = cell(length(Tamb),countTlimit);
CutLineTopXResults_vs_PSCTamb = cell(length(Tamb),countTlimit);
CutLineBackXResults_vs_PSCTamb = cell(length(Tamb),countTlimit);

% Absorption, Heat Production, Convective and radiative heat flux
AbsorptionCoeff_vs_PSCTamb = zeros(length(Tamb),countTlimit);
HeatProduction_vs_PSCTamb = zeros(length(Tamb),countTlimit);
ConvAveTop_vs_PSCTamb = zeros(length(Tamb),countTlimit);
ConvAveBack_vs_PSCTamb = zeros(length(Tamb),countTlimit);
ConvAveSide_vs_PSCTamb = zeros(length(Tamb),countTlimit);
RadAveTop_vs_PSCTamb = zeros(length(Tamb),countTlimit);
RadAveBack_vs_PSCTamb = zeros(length(Tamb),countTlimit);
RadAveSide_vs_PSCTamb = zeros(length(Tamb),countTlimit);

% Loop over environmental conditions
for idxTamb = 1:length(Tamb)
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
    EffRefAtT = effRef*(1+(n*consts.kB*EffRefT)/(VocRefG*consts.q)*log((Irradiance)/(EffRefG)))*(1+EffShiftVsTPVK*(Tcell-EffRefT)/100);
    effEref = EffRefAtT/totAbsCoeff;

    % Heat production
    QprodCoeff = (1-effEref)*totAbsCoeff;

    % Irradiance effect
    Qcell = Irradiance*QprodCoeff; % W/m^2

    % Build Comsol model
    [model,TcellAve,TcellAveTop,TcellAveBack,TmodAve,TcellMax,Conv_AveTop,Conv_AveBack,Conv_AveSide,Rad_AveTop,Rad_AveBack,Rad_AveSide,CutLineCells,CutLineTop,CutLineBack,CutLineCellsX,CutLineTopX,CutLineBackX]=PSCpanelTcomsolModel_v8(Tamb(idxTamb),v_wind,cellDistance,Qcell,pathFolder,Tcell);
    
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

    % Cell and Module Temperature
    TaveCells_vs_PSCTamb(idxTamb,countT) = Tcell;
    TAveTop_vs_PSCTamb(idxTamb,countT) = TcellAveTop;
    TAveBack_vs_PSCTamb(idxTamb,countT) = TcellAveBack;
    TaveMod_vs_PSCTamb(idxTamb,countT) = TmodAve;
    TmaxCells_vs_PSCTamb(idxTamb,countT) = TcellMax;
    % Efficiency
    Efficiency_vs_PSCTamb(idxTamb,countT)= EffRefAtT;
    % Cutlines
    CutLineCellsResults_vs_PSCTamb{idxTamb,countT} = CutLineCells;
    CutLineTopResults_vs_PSCTamb{idxTamb,countT} = CutLineTop;
    CutLineBackResults_vs_PSCTamb{idxTamb,countT} = CutLineBack;
    CutLineCellsXResults_vs_PSCTamb{idxTamb,countT} = CutLineCellsX;
    CutLineTopXResults_vs_PSCTamb{idxTamb,countT} = CutLineTopX;
    CutLineBackXResults_vs_PSCTamb{idxTamb,countT} = CutLineBackX;
    % Others
    AbsorptionCoeff_vs_PSCTamb(idxTamb,countT) = totAbsCoeff;
    HeatProduction_vs_PSCTamb(idxTamb,countT) = Qcell;
    ConvAveTop_vs_PSCTamb(idxTamb,countT) = Conv_AveTop;
    ConvAveBack_vs_PSCTamb(idxTamb,countT) = Conv_AveBack;
    ConvAveSide_vs_PSCTamb(idxTamb,countT) = Conv_AveSide;
    RadAveTop_vs_PSCTamb(idxTamb,countT) = Rad_AveTop;
    RadAveBack_vs_PSCTamb(idxTamb,countT) = Rad_AveBack;
    RadAveSide_vs_PSCTamb(idxTamb,countT) = Rad_AveSide;

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

varsRemove = {'ME','model'};
clear(varsRemove{:});

% Save workspace
save([pathFolder,'\workspace.mat'])

% Solar irradiance 

% Folder
pathFolder ='.\PSCIrr';

% Error log
FolderErrorLog = [pathFolder,'\errorLog'];
mkdir(FolderErrorLog);

% Initial temperature 
Tstc = 293.15; % K

% Weather conditions
Tamb = 20; % C
Irradiance = [100 200 300 400 500 600 700 800 900 1000 1100 1200]; % W/m^2
v_wind = 1; % m/s (unitless in comsol model)

% Panel design
cellDistance = 1; %cm

% Temperature convergence
TendCriterion = 0.1; % K
countTlimit = 10;

%
TaveCells_vs_PSCIrradiance = zeros(length(Irradiance),countTlimit);
TAveTop_vs_PSCIrradiance = zeros(length(Irradiance),countTlimit);
TAveBack_vs_PSCIrradiance = zeros(length(Irradiance),countTlimit);
TaveMod_vs_PSCIrradiance = zeros(length(Irradiance),countTlimit);
TmaxCells_vs_PSCIrradiance = zeros(length(Irradiance),countTlimit);
%
Efficiency_vs_PSCIrradiance = zeros(length(Irradiance),countTlimit);
%
CutLineCellsResults_vs_PSCIrradiance = cell(length(Irradiance),countTlimit);
CutLineTopResults_vs_PSCIrradiance = cell(length(Irradiance),countTlimit);
CutLineBackResults_vs_PSCIrradiance = cell(length(Irradiance),countTlimit);
CutLineCellsXResults_vs_PSCIrradiance = cell(length(Irradiance),countTlimit);
CutLineTopXResults_vs_PSCIrradiance = cell(length(Irradiance),countTlimit);
CutLineBackXResults_vs_PSCIrradiance = cell(length(Irradiance),countTlimit);
%
AbsorptionCoeff_vs_PSCIrradiance = zeros(length(Irradiance),countTlimit);
HeatProduction_vs_PSCIrradiance = zeros(length(Irradiance),countTlimit);

for idxIrr = 1:length(Irradiance)
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

    effEref = EffRefAtT/totAbsCoeff;

    % Heat production
    QprodCoeff = (1-effEref)*totAbsCoeff;

    % Irradiance effect
    Qcell = Irradiance(idxIrr)*QprodCoeff; % W/m^2

    % Build Comsol model
    [model,TcellAve,TcellAveTop,TcellAveBack,TmodAve,TcellMax,Conv_AveTop,Conv_AveBack,Conv_AveSide,Rad_AveTop,Rad_AveBack,Rad_AveSide,CutLineCells,CutLineTop,CutLineBack,CutLineCellsX,CutLineTopX,CutLineBackX]=PSCpanelTcomsolModel_v8(Tamb,v_wind,cellDistance,Qcell,pathFolder,Tcell);

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
    %
    TaveCells_vs_PSCIrradiance(idxIrr,countT) = Tcell;
    TaveMod_vs_PSCIrradiance(idxIrr,countT) = TmodAve;
    TAveTop_vs_PSCIrradiance(idxIrr,countT) = TcellAveTop;
    TAveBack_vs_PSCIrradiance(idxIrr,countT) = TcellAveBack;
    TmaxCells_vs_PSCIrradiance(idxIrr,countT) = TcellMax;
    %
    Efficiency_vs_PSCIrradiance(idxIrr,countT) = EffRefAtT;
    %
    CutLineCellsResults_vs_PSCIrradiance{idxIrr,countT} = CutLineCells;
    CutLineTopResults_vs_PSCIrradiance{idxIrr,countT} = CutLineTop;
    CutLineBackResults_vs_PSCIrradiance{idxIrr,countT} = CutLineBack;
    CutLineCellsXResults_vs_PSCIrradiance{idxIrr,countT} = CutLineCellsX;
    CutLineTopXResults_vs_PSCIrradiance{idxIrr,countT} = CutLineTopX;
    CutLineBackXResults_vs_PSCIrradiance{idxIrr,countT} = CutLineBackX;
    %
    AbsorptionCoeff_vs_PSCIrradiance(idxIrr,countT) = totAbsCoeff;
    HeatProduction_vs_PSCIrradiance(idxIrr,countT) = Qcell;

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

varsRemove = {'ME','model'};
clear(varsRemove{:});

% Save workspace
save([pathFolder,'\workspace.mat'])


% Wind speed %

% Folder
pathFolder ='.\PSCwind';

% Error log
FolderErrorLog = [pathFolder,'\errorLog'];
mkdir(FolderErrorLog);

% Initial temperature 
Tstc = 293.15; % K

% Weather conditions
Tamb = 20; % C
Irradiance = 800; % W/m^2
v_wind = [0 1 2 3 4 5 6 7 8 9 10]; % m/s (unitless in comsol model)

% Panel design
cellDistance = 1; %cm

% Temperature convergence
TendCriterion = 0.1; % K
countTlimit = 10;

%
TaveCells_vs_PSCWindspeed = zeros(length(v_wind),countTlimit);
TAveTop_vs_PSCWindspeed = zeros(length(v_wind),countTlimit);
TAveBack_vs_PSCWindspeed = zeros(length(v_wind),countTlimit);
TaveMod_vs_PSCWindspeed = zeros(length(v_wind),countTlimit);
TmaxCells_vs_PSCWindspeed = zeros(length(v_wind),countTlimit);
%
Efficiency_vs_PSCWindspeed = zeros(length(v_wind), countTlimit);
%
CutLineCellsResults_vs_PSCWindspeed = cell(length(v_wind),countTlimit);
CutLineTopResults_vs_PSCWindspeed = cell(length(v_wind),countTlimit);
CutLineBackResults_vs_PSCWindspeed = cell(length(v_wind),countTlimit);
CutLineCellsXResults_vs_PSCWindspeed = cell(length(v_wind),countTlimit);
CutLineTopXResults_vs_PSCWindspeed = cell(length(v_wind),countTlimit);
CutLineBackXResults_vs_PSCWindspeed = cell(length(v_wind),countTlimit);
%
AbsorptionCoeff_vs_PSCWindspeed = zeros(length(v_wind),countTlimit);
HeatProduction_vs_PSCWindspeed = zeros(length(v_wind),countTlimit);

for idxv_wind = 1:length(v_wind)
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

    % Electrical efficience
    EffRefAtT = effRef*(1+(n*consts.kB*EffRefT)/(VocRefG*consts.q)*log((Irradiance)/(EffRefG)))*(1+EffShiftVsTPVK*(Tcell-EffRefT)/100);

    effEref = EffRefAtT/totAbsCoeff;

    % Heat production
    QprodCoeff = (1-effEref)*totAbsCoeff;

    % Irradiance effect
    Qcell = Irradiance*QprodCoeff; % W/m^2

    % Build Comsol model
    [model,TcellAve,TcellAveTop,TcellAveBack,TmodAve,TcellMax,Conv_AveTop,Conv_AveBack,Conv_AveSide,Rad_AveTop,Rad_AveBack,Rad_AveSide,CutLineCells,CutLineTop,CutLineBack,CutLineCellsX,CutLineTopX,CutLineBackX]=PSCpanelTcomsolModel_v8(Tamb,v_wind(idxv_wind),cellDistance,Qcell,pathFolder,Tcell);

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
    %
    TaveCells_vs_PSCWindspeed(idxv_wind,countT) = Tcell;
    TAveTop_vs_PSCWindspeed(idxv_wind,countT) = TcellAveTop;
    TAveBack_vs_PSCWindspeed(idxv_wind,countT) = TcellAveBack;
    TaveMod_vs_PSCWindspeed(idxv_wind,countT) = TmodAve;
    TmaxCells_vs_PSCWindspeed(idxv_wind,countT) = TcellMax;
    %
    Efficiency_vs_PSCWindspeed(idxv_wind,countT) = EffRefAtT;
    %
    CutLineCellsResults_vs_PSCWindspeed{idxv_wind,countT} = CutLineCells;
    CutLineTopResults_vs_PSCWindspeed{idxv_wind,countT} = CutLineTop;
    CutLineBackResults_vs_PSCWindspeed{idxv_wind,countT} = CutLineBack;
    CutLineCellsXResults_vs_PSCWindspeed{idxv_wind,countT} = CutLineCellsX;
    CutLineTopXResults_vs_PSCWindspeed{idxv_wind,countT} = CutLineTopX;
    CutLineBackXResults_vs_PSCWindspeed{idxv_wind,countT} = CutLineBackX;
    %
    AbsorptionCoeff_vs_PSCWindspeed(idxv_wind,countT) = totAbsCoeff;
    HeatProduction_vs_PSCWindspeed(idxv_wind,countT) = Qcell;

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

varsRemove = {'ME','model'};
clear(varsRemove{:});

% Save workspace
save([pathFolder,'\workspace.mat'])

% SILICON %
errorCountLimit=5;
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

% Total efficiencies
EffRefT = 298.15; % K
EffRefG = 1000; % W/m^2
VocRefG = 0.6; % W
n=1.3;
% Si
effRef = 0.2;
EffShiftVsTSi = -0.35; % %/K

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

% Ambient temperature % 

% Folder
pathFolder ='.\SiTamb';
% 
% Error log
FolderErrorLog = [pathFolder,'\errorLog'];
mkdir(FolderErrorLog);

% Initial temperature 
Tstc = 293.15; % K

% Weather conditions
Tamb = 20; %[-20 -10 0 10 20 30 40 50]; % C
Irradiance = 800; % W/m^2
v_wind = 1; % m/s (unitless in comsol model)

% Panel design
cellDistance = 1; %cm

% Temperature convergence
TendCriterion = 0.1; % K
countTlimit = 10;

TaveCells_vs_SiTamb = zeros(length(Tamb),countTlimit);
TAveTop_vs_SiTamb = zeros(length(Tamb),countTlimit);
TAveBack_vs_SiTamb = zeros(length(Tamb),countTlimit);
TaveMod_vs_SiTamb = zeros(length(Tamb),countTlimit);
TmaxCells_vs_SiTamb = zeros(length(Tamb),countTlimit);
%
Efficiency_vs_SiTamb = zeros(length(Tamb),countTlimit);
CutLineCellsResults_vs_SiTamb = cell(length(Tamb),countTlimit);
CutLineTopResults_vs_SiTamb = cell(length(Tamb),countTlimit);
CutLineBackResults_vs_SiTamb = cell(length(Tamb),countTlimit);
CutLineCellsXResults_vs_SiTamb = cell(length(Tamb),countTlimit);
CutLineTopXResults_vs_SiTamb = cell(length(Tamb),countTlimit);
CutLineBackXResults_vs_SiTamb = cell(length(Tamb),countTlimit);
%
AbsorptionCoeff_vs_SiTamb = zeros(length(Tamb),countTlimit);
HeatProduction_vs_SiTamb = zeros(length(Tamb),countTlimit);
ConvAveTop_vs_SiTamb = zeros(length(Tamb),countTlimit);
ConvAveBack_vs_SiTamb = zeros(length(Tamb),countTlimit);
ConvAveSide_vs_SiTamb = zeros(length(Tamb),countTlimit);
RadAveTop_vs_SiTamb = zeros(length(Tamb),countTlimit);
RadAveBack_vs_SiTamb = zeros(length(Tamb),countTlimit);
RadAveSide_vs_SiTamb = zeros(length(Tamb),countTlimit);

% Loop over environmental conditions
for idxTamb = 1:length(Tamb)

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
    EffRefAtT = effRef*(1+(n*consts.kB*EffRefT)/(VocRefG*consts.q)*log((Irradiance)/(EffRefG)))*(1+EffShiftVsTSi*(Tcell-EffRefT)/100);

    effEref = EffRefAtT/totAbsCoeff;

    % Heat production
    QprodCoeff = (1-effEref)*totAbsCoeff;

    % Irradiance effect
    Qcell = Irradiance*QprodCoeff; % W/m^2

    % Build Comsol model
    [model,TcellAve,TcellAveTop,TcellAveBack,TmodAve,TcellMax,Conv_AveTop,Conv_AveBack,Conv_AveSide,Rad_AveTop,Rad_AveBack,Rad_AveSide,CutLineCells,CutLineTop,CutLineBack,CutLineCellsX,CutLineTopX,CutLineBackX] = SipanelTcomsolModel_v6(Tamb(idxTamb),v_wind,cellDistance,Qcell,pathFolder,Tcell);
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
    TaveCells_vs_SiTamb(idxTamb,countT) = Tcell;
    TAveTop_vs_SiTamb(idxTamb,countT) = TcellAveTop;
    TAveBack_vs_SiTamb(idxTamb,countT) = TcellAveBack;
    TaveMod_vs_SiTamb(idxTamb,countT) = TmodAve;
    TmaxCells_vs_SiTamb(idxTamb,countT) = TcellMax;
    %
    Efficiency_vs_SiTamb(idxTamb,countT) = EffRefAtT;
    %
    CutLineCellsResults_vs_SiTamb{idxTamb,countT} = CutLineCells;
    CutLineTopResults_vs_SiTamb{idxTamb,countT} = CutLineTop;
    CutLineBackResults_vs_SiTamb{idxTamb,countT} = CutLineBack;
    CutLineCellsXResults_vs_SiTamb{idxTamb,countT} = CutLineCellsX;
    CutLineTopXResults_vs_SiTamb{idxTamb,countT} = CutLineTopX;
    CutLineBackXResults_vs_SiTamb{idxTamb,countT} = CutLineBackX;
    %
    AbsorptionCoeff_vs_SiTamb(idxTamb,countT)=totAbsCoeff;
    HeatProduction_vs_SiTamb(idxTamb,countT) = Qcell;
    ConvAveTop_vs_SiTamb(idxTamb,countT) = Conv_AveTop;
    ConvAveBack_vs_SiTamb(idxTamb,countT) = Conv_AveBack;
    ConvAveSide_vs_SiTamb(idxTamb,countT) = Conv_AveSide;
    RadAveTop_vs_SiTamb(idxTamb,countT) = Rad_AveTop;
    RadAveBack_vs_SiTamb(idxTamb,countT) = Rad_AveBack;
    RadAveSide_vs_SiTamb(idxTamb,countT) = Rad_AveSide;


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

varsRemove = {'ME','model'};
clear(varsRemove{:});

% Save workspace
save([pathFolder,'\workspace.mat'])


% % Solar irradiance 

% % Folder
pathFolder ='.\SiIrr';

% Error log
FolderErrorLog = [pathFolder,'\errorLog'];
mkdir(FolderErrorLog);

% Initial temperature 
Tstc = 293.15; % K

% Weather conditions
Tamb = 20; % C
Irradiance = [100 200 300 400 500 600 700 800 900 1000 1100 1200]; % W/m^2
v_wind = 1; % m/s (unitless in comsol model)

% Panel design
cellDistance = 1; %cm
% Temperature convergence
TendCriterion = 0.1; % K
countTlimit = 10;

%
TaveCells_vs_SiIrradiance = zeros(length(Irradiance),countTlimit);
TAveTop_vs_SiIrradiance = zeros(length(Irradiance),countTlimit);
TAveBack_vs_SiAll = zeros(length(Irradiance),countTlimit);
TaveMod_vs_SiIrradiance = zeros(length(Irradiance),countTlimit);
TmaxCells_vs_SiIrradiance = zeros(length(Irradiance),countTlimit);
%
Efficiency_vs_SiIrradiance = zeros(length(Irradiance),countTlimit);
%
CutLineCellsResults_vs_SiIrradiance = cell(length(Irradiance),countTlimit);
CutLineTopResults_vs_SiIrradiance = cell(length(Irradiance),countTlimit);
CutLineBackResults_vs_SiIrradiance = cell(length(Irradiance),countTlimit);
CutLineCellsXResults_vs_SiIrradiance = cell(length(Irradiance),countTlimit);
CutLineTopXResults_vs_SiIrradiance = cell(length(Irradiance),countTlimit);
CutLineBackXResults_vs_SiIrradiance = cell(length(Irradiance),countTlimit);
%
AbsorptionCoeff_vs_SiIrradiance = cell(length(Irradiance),countTlimit);
HeatProduction_vs_SiIrradiance = cell(length(Irradiance),countTlimit);

for idxIrr = 1:length(Irradiance)
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
    EffRefAtT = effRef*(1+(n*consts.kB*EffRefT)/(VocRefG*consts.q)*log((Irradiance(idxIrr))/(EffRefG)))*(1+EffShiftVsTSi*(Tcell-EffRefT)/100);

    effEref = EffRefAtT/totAbsCoeff;

    % Heat production
    QprodCoeff = (1-effEref)*totAbsCoeff;

    % Irradiance effect
    Qcell = Irradiance(idxIrr)*QprodCoeff; % W/m^2

    % Build Comsol model
    [model,TcellAve,TcellAveTop,TcellAveBack,TmodAve,TcellMax,Conv_AveTop,Conv_AveBack,Conv_AveSide,Rad_AveTop,Rad_AveBack,Rad_AveSide,CutLineCells,CutLineTop,CutLineBack,CutLineCellsX,CutLineTopX,CutLineBackX] = SipanelTcomsolModel_v6(Tamb,v_wind,cellDistance,Qcell,pathFolder,Tcell);
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
    TaveCells_vs_SiIrradiance(idxIrr,countT) = Tcell;
    TAveTop_vs_SiIrradiance(idxIrr,countT) = TcellAveTop;
    TAveBack_vs_SiAll(idxIrr,countT) = TcellAveBack;
    TaveMod_vs_SiIrradiance(idxIrr,countT) = TmodAve;
    TmaxCells_vs_SiIrradiance(idxIrr,countT) = TcellMax;
    %
    Efficiency_vs_SiIrradiance(idxIrr,countT) = EffRefAtT;
    %
    CutLineCellsResults_vs_SiIrradiance{idxIrr,countT} = CutLineCells;
    CutLineTopResults_vs_SiIrradiance{idxIrr,countT} = CutLineTop;
    CutLineBackResults_vs_SiIrradiance{idxIrr,countT} = CutLineBack;
    CutLineCellsXResults_vs_SiIrradiance{idxIrr,countT} = CutLineCellsX;
    CutLineTopXResults_vs_SiIrradiance{idxIrr,countT} = CutLineTopX;
    CutLineBackXResults_vs_SiIrradiance{idxIrr,countT} = CutLineBackX;
    %
    AbsorptionCoeff_vs_SiIrradiance(idxIrr,countT)=totAbsCoeff;
    HeatProduction_vs_SiIrradiance(idxIrr,countT) = Qcell;


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

varsRemove = {'ME','model'};
clear(varsRemove{:});

% Save workspace
save([pathFolder,'\workspace.mat'])


% Wind speed

% Folder
pathFolder ='.\Siwind';

% Error log
FolderErrorLog = [pathFolder,'\errorLog'];
mkdir(FolderErrorLog);

% Initial temperature 
Tstc = 293.15; % K

% Weather conditions
Tamb = 20; % C
Irradiance = 800; % W/m^2
v_wind = [0 1 2 3 4 5 6 7 8 9 10]; % m/s (unitless in comsol model)

% Panel design
cellDistance = 1; %cm
% Temperature convergence
TendCriterion = 0.1; % K
countTlimit = 10;

%
TaveCells_vs_SiWindspeed = zeros(length(v_wind),countTlimit);
TAveTop_vs_SiWindspeed = zeros(length(v_wind),countTlimit);
TAveBack_vs_SiWindspeed = zeros(length(v_wind),countTlimit);
TaveMod_vs_SiWindspeed = zeros(length(v_wind),countTlimit);
TmaxCells_vs_SiWindspeed = zeros(length(v_wind),countTlimit);
%
Efficiency_vs_SiWindspeed = zeros(length(v_wind),countTlimit);
%
CutLineCellsResults_vs_SiWindspeed = cell(length(v_wind),countTlimit);
CutLineTopResults_vs_SiWindspeed = cell(length(v_wind),countTlimit);
CutLineBackResults_vs_SiWindspeed = cell(length(v_wind),countTlimit);
CutLineCellsXResults_vs_SiWindspeed = cell(length(v_wind),countTlimit);
CutLineTopXResults_vs_SiWindspeed = cell(length(v_wind),countTlimit);
CutLineBackXResults_vs_SiWindspeed = cell(length(v_wind),countTlimit);
%
AbsorptionCoeff_vs_SiWindspeed = cell(length(v_wind),countTlimit);
HeatProduction_vs_SiWindspeed = cell(length(v_wind),countTlimit);


for idxv_wind = 1:length(v_wind)
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
    EffRefAtT = effRef*(1+(n*consts.kB*EffRefT)/(VocRefG*consts.q)*log((Irradiance)/(EffRefG)))*(1+EffShiftVsTSi*(Tcell-EffRefT)/100);

    effEref = EffRefAtT/totAbsCoeff;

    % Heat production
    QprodCoeff = (1-effEref)*totAbsCoeff;

    % Irradiance effect
    Qcell = Irradiance*QprodCoeff; % W/m^2

    % Build Comsol model
    [model,TcellAve,TcellAveTop,TcellAveBack,TmodAve,TcellMax,Conv_AveTop,Conv_AveBack,Conv_AveSide,Rad_AveTop,Rad_AveBack,Rad_AveSide,CutLineCells,CutLineTop,CutLineBack,CutLineCellsX,CutLineTopX,CutLineBackX] = SipanelTcomsolModel_v6(Tamb,v_wind(idxv_wind),cellDistance,Qcell,pathFolder,Tcell);
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
    %
    TaveCells_vs_SiWindspeed(idxv_wind,countT) = Tcell;
    TAveTop_vs_SiWindspeed(idxv_wind,countT) = TcellAveTop;
    TAveBack_vs_SiAll(idxv_wind,countT) = TcellAveBack;
    TaveMod_vs_SiWindspeed(idxv_wind,countT) = TmodAve;
    TmaxCells_vs_SiWindspeed(idxv_wind,countT) = TcellMax;
    %
    Efficiency_vs_SiWindspeed(idxv_wind,countT) = EffRefAtT;
    %
    CutLineCellsResults_vs_SiWindspeed{idxv_wind,countT} = CutLineCells;
    CutLineTopResults_vs_SiWindspeed{idxv_wind,countT} = CutLineTop;
    CutLineBackResults_vs_SiWindspeed{idxv_wind,countT} = CutLineBack;
    CutLineCellsXResults_vs_SiWindspeed{idxv_wind,countT} = CutLineCellsX;
    CutLineTopXResults_vs_SiWindspeed{idxv_wind,countT} = CutLineTopX;
    CutLineBackXResults_vs_SiWindspeed{idxv_wind,countT} = CutLineBackX;
    %
    AbsorptionCoeff_vs_SiWindspeed(idxv_wind,countT)= totAbsCoeff;
    HeatProduction_vs_SiWindspeed(idxv_wind,countT) = Qcell;

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

varsRemove = {'ME','model'};
clear(varsRemove{:});

% Save workspace
save([pathFolder,'\workspace.mat'])

