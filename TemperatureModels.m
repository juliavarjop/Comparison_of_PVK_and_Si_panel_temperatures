% This code is used to determine the parameters for Sandia, Faiman, Mattei, TRNSYS, and PVsyst temperature models

% Use SaveResults.m to store the simulated data required for determining the temperature model parameters

data = TaveCellsvsSiIrrWind;
data2 = AbsCoeffvsSiIrrWind; 

% Model inputs
Tcell = data(:,2:7);
totAbsCoeff = data2(:,2:7);
WS = data(:,1);
Tamb = 20;
E = 200:200:1200;
at=totAbsCoeff;
n=0.2;
beta=-0.0035;
Tr=25;

% Sandia Model Parameters

y=log((Tcell-Tamb)./E);

y1=y(:,1);
y2=y(:,2);
y3=y(:,3);
y4=y(:,4);
y5=y(:,5);
y6=y(:,6);

Xdata_WS = repmat(WS,size(E));

Xdata_WS = repmat(WS,size(transpose(E)));

x1 = Xdata_WS>=0;
x2 = Xdata_WS<=4;
x3 = Xdata_WS>=4;
x4 = Xdata_WS<=10;

TableYdata=y; 

Ydata = TableYdata; 

Ydata = TableYdata(:); 
% 
Xdata_WS = repmat(WS,size(E)); 

Xdata_WS = repmat(WS,size(transpose(E))); 

% fit between 0-4 m/s
Sandia_p_0to4 = polyfit(Xdata_WS(x1 & x2),Ydata(x1 & x2),1);

% fit between 4-10 m/s 
Sandia_p_4to10 = polyfit(Xdata_WS(x3 & x4),Ydata(x3 & x4),1);

% fit between 0-10 m/s
Sandia_p_0to10 = polyfit(Xdata_WS(x1 & x4),Ydata(x1 & x4),1)

a=Sandia_p_0to10(:,2);
b=Sandia_p_0to10(:,1);

% Faiman Model Parameters

y = E./(Tcell-Tamb);

y1=y(:,1);
y2=y(:,2);
y3=y(:,3);
y4=y(:,4);
y5=y(:,5);
y6=y(:,6);

TableYdata=y; 

Ydata = TableYdata; 

Ydata = TableYdata(:); 

Xdata_WS = repmat(WS,size(E)); 

Xdata_WS = repmat(WS,size(transpose(E))); 

x1 = Xdata_WS>=0;
x4 = Xdata_WS<=8;

% fit between 0-10 m/s
Faiman_p_0to10 = polyfit(Xdata_WS(x1 & x4),Ydata(x1 & x4),1)

U_L1=Faiman_p_0to10(:,1);
U_L0=Faiman_p_0to10(:,2);

% Mattei Model Parameters

y=(E.*(at-n-beta*n.*(Tcell-Tr)))./(Tcell-Tamb);

% Cell Distance 1 cm
y1=y(:,1);
y2=y(:,2);
y3=y(:,3);
y4=y(:,4);
y5=y(:,5);
y6=y(:,6);

TableYdata=y; 

Ydata = TableYdata; 

Ydata = TableYdata(:); 

Xdata_WS = repmat(WS,size(E)); 

Xdata_WS = repmat(WS,size(transpose(E)));

x1 = Xdata_WS>=0;
x4 = Xdata_WS<=10;

% fit between 0-10 m/s
Mattei_p_0to10 = polyfit(Xdata_WS(x1 & x4),Ydata(x1 & x4),1)

B=Mattei_p_0to10(:,1);
A=Mattei_p_0to10(:,2);
 
% TRNSYS Parameters

y=(at.*E.*(1-(n./at)))./(Tcell-Tamb);  

y1=y(:,1);
y2=y(:,2);
y3=y(:,3);
y4=y(:,4);
y5=y(:,5);
y6=y(:,6);

TableYdata=y; 

Ydata = TableYdata; 

Ydata = TableYdata(:); 

Xdata_WS = repmat(WS,size(E)); 

Xdata_WS = repmat(WS,size(transpose(E))); 

x1 = Xdata_WS>=0;
x4 = Xdata_WS<=8;

% fit between 0-10 m/s
TRNSYS_p_0to10 = polyfit(Xdata_WS(x1 & x4),Ydata(x1 & x4),1)

U_Loss1=TRNSYS_p_0to10(:,1);
U_Loss2=TRNSYS_p_0to10(:,2);

U_Loss=U_Loss1+U_Loss2;

% PVsyst Model Parameters
y = (E.*at.*(1-n))./(Tcell-Tamb);

% Cell Distance 1 cm
y1=y(:,1);
y2=y(:,2);
y3=y(:,3);
y4=y(:,4);
y5=y(:,5);
y6=y(:,6);

TableYdata=y; 

Ydata = TableYdata; 

Ydata = TableYdata(:); 

Xdata_WS = repmat(WS,size(E)); 

Xdata_WS = repmat(WS,size(transpose(E)));

x1 = Xdata_WS>=0;
x4 = Xdata_WS<=10;

% fit between 0-10 m/s
PVsyst_p_0to10 = polyfit(Xdata_WS(x1 & x4),Ydata(x1 & x4),1)

U1=PVsyst_p_0to10(:,1);
U0=PVsyst_p_0to10(:,2);

