
% The data simulated in varying combinations of ambient conditions (main_PVKpanelT_simpleComparison_AllConditions.m)
% load('ExampleDataSilicon.mat')

indDataPoint = 1;
indCellDistance = 1;
indkcoeff=3; % The constant used in the convective heat transfer coefficient equation h = kcoeff * v_wind + 2.8 (here kcoeff is 3 Ws/m^3K)
for indWS = 1:length(v_wind) 
for indIrr = 1:length(Irradiance)
for indTamb = 1:length(Tamb)

    ind = find(TaveCells_vs_SiAll(indTamb,indIrr,indWS,indCellDistance,indkcoeff,:),1,'last');
    ind1 = find(AbsorptionCoeff_vs_SiAll(indTamb,indIrr,indWS,indCellDistance,indkcoeff,:),1,'last');

    indT = find(~cellfun('isempty', CutLineBackResults_vs_SiAll(indTamb,indIrr,indWS,indCellDistance,indkcoeff,:)), 1, 'last'); 
    indx = find(~cellfun('isempty', CutLineBackXResults_vs_SiAll(indTamb,indIrr,indWS,indCellDistance,indkcoeff,:)), 1, 'last'); 
    
    x=CutLineBackXResults_vs_SiAll{indTamb,indIrr,indWS,indCellDistance,indkcoeff,indx}*100;
    y=CutLineBackResults_vs_SiAll{indTamb,indIrr,indWS,indCellDistance,indkcoeff,indT};
   
    n = 1;
    m = 1;
    limitMIN = 11;
    limitMAX = 12;
    listMIN = [];
    listMAX = [];

    for i = 1:limitMIN
    idx = find(x>=n*(1+15.6) & x<=(n+1)*1+n*15.6);
    yidx = y(idx);
    n = n + 1;

        listMIN = [listMIN;{yidx}];
        if n > limitMIN
        end

    end

    for j = 1:limitMAX
    idx2 = find(x>=m*1+(m-1)*15.6 & x<=m*(1+15.6));
    yidx = y(idx2);
    m = m + 1;

        listMAX = [listMAX;{yidx}];
        if m > limitMAX
        end
    end
    
    % Store results 
    
    % The maximum temperature of the module back surface in varying combinations of ambient conditions
    TmaxBackvsSiAllCond(indDataPoint,1) = Tamb(indTamb);
    TmaxBackvsSiAllCond(indDataPoint,2) = Irradiance(indIrr);
    TmaxBackvsSiAllCond(indDataPoint,3) = v_wind(indWS);
    TmaxBackvsSiAllCond(indDataPoint,4)=max([listMAX{:}])-273.15;

    % The minimum temperature of the module back surface in varying combinations of ambient conditions
    TminBackvsSiAllCond(indDataPoint,1) = Tamb(indTamb);
    TminBackvsSiAllCond(indDataPoint,2) = Irradiance(indIrr);
    TminBackvsSiAllCond(indDataPoint,3) = v_wind(indWS);
    TminBackvsSiAllCond(indDataPoint,4)=min([listMIN{:}])-273.15;
    
    % The average temperature of the cells in varying combinations of ambient conditions
    TaveCellsvsSiAllCond(indDataPoint,1) = Tamb(indTamb);
    TaveCellsvsSiAllCond(indDataPoint,2) = Irradiance(indIrr);
    TaveCellsvsSiAllCond(indDataPoint,3) = v_wind(indWS);
    TaveCellsvsSiAllCond(indDataPoint,4) = TaveCells_vs_SiAll(indTamb,indIrr,indWS,indCellDistance,indkcoeff,ind)-273.15;
    
    % The absorption coefficient in varying combinations of ambient conditions
    AbsCoeffvsSiAllCond(indDataPoint,1) = Tamb(indTamb);
    AbsCoeffvsSiAllCond(indDataPoint,2) = Irradiance(indIrr);
    AbsCoeffvsSiAllCond(indDataPoint,3) = v_wind(indWS);
    AbsCoeffvsSiAllCond(indDataPoint,4) = AbsorptionCoeff_vs_SiAll(indTamb,indIrr,indWS,indCellDistance,indkcoeff,ind);

    indDataPoint = indDataPoint+1;

end
end
end

% The maximum temperature of module back surface in varying combinations of wind speed and irradiance

filteredData1 = TmaxBackvsSiAllCond(TmaxBackvsSiAllCond(:,1) == 20, :);

TmaxBackvsSiIrrWind = nan(length(v_wind), length(Irradiance)+1);

TmaxBackvsSiIrrWind(:,1) = v_wind;

for i = 1:size(filteredData1, 1)
    idxWind = find(v_wind == filteredData1(i,3)); 
    idxIrr = find(Irradiance == filteredData1(i,2)); 
    TmaxBackvsSiIrrWind(idxWind, idxIrr+1) = filteredData1(i,4);
end

% The minimum temperature of module back surface in varying combinations of wind speed and irradiance
filteredData = TminBackvsSiAllCond(TminBackvsSiAllCond(:,1) == 20, :);

TminBackvsSiIrrWind = nan(length(v_wind), length(Irradiance)+1);

TminBackvsSiIrrWind(:,1) = v_wind;

for i = 1:size(filteredData, 1)
    idxWind = find(v_wind == filteredData(i,3)); 
    idxIrr = find(Irradiance == filteredData(i,2)); 
    TminBackvsSiIrrWind(idxWind, idxIrr+1) = filteredData(i,4);
end

% The average temperature of the cells in varying combinations of wind speed and irradiance

filteredData = TaveCellsvsSiAllCond(TaveCellsvsSiAllCond(:,1) == 20, :);

TaveCellsvsSiIrrWind = nan(length(v_wind), length(Irradiance)+1);

TaveCellsvsSiIrrWind(:,1) = v_wind;

for i = 1:size(filteredData, 1)
    idxWind = find(v_wind == filteredData(i,3)); 
    idxIrr = find(Irradiance == filteredData(i,2)); 
    TaveCellsvsSiIrrWind(idxWind, idxIrr+1) = filteredData(i,4);
end

% Absorption coefficient in varying combinations of wind speed and irradiance

filteredData = AbsCoeffvsSiAllCond(AbsCoeffvsSiAllCond(:,1) == 20, :);

AbsCoeffvsSiIrrWind = nan(length(v_wind), length(Irradiance)+1);

AbsCoeffvsSiIrrWind(:,1) = v_wind;

for i = 1:size(filteredData, 1)
    idxWind = find(v_wind == filteredData(i,3)); 
    idxIrr = find(Irradiance == filteredData(i,2)); 
    AbsCoeffvsSiIrrWind(idxWind, idxIrr+1) = filteredData(i,4);
end

% load('ExampleDataPerovskite.mat')

indDataPoint = 1;
indCellDistance = 1;
indkcoeff=3; % The constant used in the convective heat transfer coefficient equation h = kcoeff * v_wind + 2.8 (here kcoeff is 3 Ws/m^3K)
for indWS = 1:length(v_wind) 
for indIrr = 1:length(Irradiance)
for indTamb = 1:length(Tamb)

    ind = find(TaveCells_vs_PSCAll(indTamb,indIrr,indWS,indCellDistance,indkcoeff,:),1,'last');
    ind1 = find(AbsorptionCoeff_vs_PSCAll(indTamb,indIrr,indWS,indCellDistance,indkcoeff,:),1,'last');

    indT = find(~cellfun('isempty', CutLineBackResults_vs_PSCAll(indTamb,indIrr,indWS,indCellDistance,indkcoeff,:)), 1, 'last'); 
    indx = find(~cellfun('isempty', CutLineBackXResults_vs_PSCAll(indTamb,indIrr,indWS,indCellDistance,indkcoeff,:)), 1, 'last'); 
    
    x=CutLineBackXResults_vs_PSCAll{indTamb,indIrr,indWS,indCellDistance,indkcoeff,indx}*100;
    y=CutLineBackResults_vs_PSCAll{indTamb,indIrr,indWS,indCellDistance,indkcoeff,indT};
   
    n = 1;
    m = 1;
    limitMIN = 11;
    limitMAX = 12;
    listMIN = [];
    listMAX = [];

    for i = 1:limitMIN
    idx = find(x>=n*(1+15.6) & x<=(n+1)*1+n*15.6);
    yidx = y(idx);
    n = n + 1;

        listMIN = [listMIN;{yidx}];
        if n > limitMIN
        end

    end

    for j = 1:limitMAX
    idx2 = find(x>=m*1+(m-1)*15.6 & x<=m*(1+15.6));
    yidx = y(idx2);
    m = m + 1;

        listMAX = [listMAX;{yidx}];
        if m > limitMAX
        end
    end
    
    % Store results 
    
    % The maximum temperature of the module back surface in varying combinations of ambient conditions
    TmaxBackvsPSCAllCond(indDataPoint,1) = Tamb(indTamb);
    TmaxBackvsPSCAllCond(indDataPoint,2) = Irradiance(indIrr);
    TmaxBackvsPSCAllCond(indDataPoint,3) = v_wind(indWS);
    TmaxBackvsPSCAllCond(indDataPoint,4)=max([listMAX{:}])-273.15;

    % The minimum temperature of the module back surface in varying combinations of ambient conditions
    TminBackvsPSCAllCond(indDataPoint,1) = Tamb(indTamb);
    TminBackvsPSCAllCond(indDataPoint,2) = Irradiance(indIrr);
    TminBackvsPSCAllCond(indDataPoint,3) = v_wind(indWS);
    TminBackvsPSCAllCond(indDataPoint,4)=min([listMIN{:}])-273.15;
    
    % The average temperature of the cells in varying combinations of ambient conditions
    TaveCellsvsPSCAllCond(indDataPoint,1) = Tamb(indTamb);
    TaveCellsvsPSCAllCond(indDataPoint,2) = Irradiance(indIrr);
    TaveCellsvsPSCAllCond(indDataPoint,3) = v_wind(indWS);
    TaveCellsvsPSCAllCond(indDataPoint,4) = TaveCells_vs_PSCAll(indTamb,indIrr,indWS,indCellDistance,indkcoeff,ind)-273.15;
    
    % The absorption coefficient in varying combinations of ambient conditions
    AbsCoeffvsPSCAllCond(indDataPoint,1) = Tamb(indTamb);
    AbsCoeffvsPSCAllCond(indDataPoint,2) = Irradiance(indIrr);
    AbsCoeffvsPSCAllCond(indDataPoint,3) = v_wind(indWS);
    AbsCoeffvsPSCAllCond(indDataPoint,4) = AbsorptionCoeff_vs_PSCAll(indTamb,indIrr,indWS,indCellDistance,indkcoeff,ind);

    indDataPoint = indDataPoint+1;

end
end
end

% The maximum temperature of module back surface in varying combinations of wind speed and irradiance

filteredData1 = TmaxBackvsPSCAllCond(TmaxBackvsPSCAllCond(:,1) == 20, :);

TmaxBackvsPSCIrrWind = nan(length(v_wind), length(Irradiance)+1);

TmaxBackvsPSCIrrWind(:,1) = v_wind;

for i = 1:size(filteredData1, 1)
    idxWind = find(v_wind == filteredData1(i,3)); 
    idxIrr = find(Irradiance == filteredData1(i,2)); 
    TmaxBackvsPSCIrrWind(idxWind, idxIrr+1) = filteredData1(i,4);
end

% The minimum temperature of module back surface in varying combinations of wind speed and irradiance
filteredData = TminBackvsPSCAllCond(TminBackvsPSCAllCond(:,1) == 20, :);

TminBackvsPSCIrrWind = nan(length(v_wind), length(Irradiance)+1);

TminBackvsPSCIrrWind(:,1) = v_wind;

for i = 1:size(filteredData, 1)
    idxWind = find(v_wind == filteredData(i,3)); 
    idxIrr = find(Irradiance == filteredData(i,2)); 
    TminBackvsPSCIrrWind(idxWind, idxIrr+1) = filteredData(i,4);
end

% The average temperature of the cells in varying combinations of wind speed and irradiance

filteredData = TaveCellsvsPSCAllCond(TaveCellsvsPSCAllCond(:,1) == 20, :);

TaveCellsvsPSCIrrWind = nan(length(v_wind), length(Irradiance)+1);

TaveCellsvsPSCIrrWind(:,1) = v_wind;

for i = 1:size(filteredData, 1)
    idxWind = find(v_wind == filteredData(i,3)); 
    idxIrr = find(Irradiance == filteredData(i,2)); 
    TaveCellsvsPSCIrrWind(idxWind, idxIrr+1) = filteredData(i,4);
end

% Absorption coefficient in varying combinations of wind speed and irradiance

filteredData = AbsCoeffvsPSCAllCond(AbsCoeffvsPSCAllCond(:,1) == 20, :);

AbsCoeffvsPSCIrrWind = nan(length(v_wind), length(Irradiance)+1);

AbsCoeffvsPSCIrrWind(:,1) = v_wind;

for i = 1:size(filteredData, 1)
    idxWind = find(v_wind == filteredData(i,3)); 
    idxIrr = find(Irradiance == filteredData(i,2)); 
    AbsCoeffvsPSCIrrWind(idxWind, idxIrr+1) = filteredData(i,4);
end
