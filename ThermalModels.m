% Use SaveResults.m to create the required files

data = TminBackvsPSCIrrWind; 
data2 = AbsCoeffvsPSCIrrWind; 

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

% Sandia Model

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
p_0to4 = polyfit(Xdata_WS(x1 & x2),Ydata(x1 & x2),1)

% fit between 4-10 m/s 
p_4to10 = polyfit(Xdata_WS(x3 & x4),Ydata(x3 & x4),1)

% fit between 0-10 m/s
p_0to10 = polyfit(Xdata_WS(x1 & x4),Ydata(x1 & x4),1)

a=p_0to10(:,2);
b=p_0to10(:,1);

% Faiman Model

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
p_0to10 = polyfit(Xdata_WS(x1 & x4),Ydata(x1 & x4),1)

U_L1=p_0to10(:,1);
U_L0=p_0to10(:,2);

% Mattei Model

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
p_0to10 = polyfit(Xdata_WS(x1 & x4),Ydata(x1 & x4),1)

B=p_0to10(:,1);
A=p_0to10(:,2);

% Skoplaki

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
p_0to10 = polyfit(Xdata_WS(x1 & x4),Ydata(x1 & x4),1)

U_Loss1=p_0to10(:,1);
U_Loss2=p_0to10(:,2);

U_Loss=U_Loss1+U_Loss2;


% PVsyst Model
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
p_0to10 = polyfit(Xdata_WS(x1 & x4),Ydata(x1 & x4),1)

U1=p_0to10(:,1);
U0=p_0to10(:,2);

%% FIGURES

% % plot colors for perovskite
redShades = [0.4, 0, 0;
            0.6, 0, 0; 
            0.8, 0, 0; 
            1, 0, 0;
            1, 0.2, 0.2; 
            1, 0.4, 0.4; 
            1, 0.6, 0.6;];

% % plot colors for silicon
blueShades = [0 0.1 0.3;
              0 0.1 0.5;
              0 0.2 0.7;
              0 0.3 0.7;
              0.1 0.3 0.9;
              0.1 0.5 0.9
              0.2 0.5 1;];

width = 625; 
height = 500;  

figure('Position', [100, 100, width, height]);

% Sandia parameters for silicon panel

% % Line shows between 0-4 m/s
% WS_0to4 = WS(WS>=0 & WS<=4);
% y_0to4 = p_0to4(:,1) * WS_0to4+ p_0to4(:,2);
% 
% % Line shows between 4-10 m/s
% WS_4to10 = WS(WS>=4 & WS<=10);
% y_4to10 = p_4to10(:,1) * WS_4to10+ p_4to10(:,2);

% Line shows between 0-10 m/s
WS_0to10 = WS(WS>=0 & WS<=10);
y_0to10 = p_0to10(:,1) * WS_0to10+ p_0to10(:,2);

% Silicon panel
% plot1=plot(WS,y1,'o','Color',blueShades(6,:),'MarkerSize',7,'LineWidth',1.5);
% hold on
% plot2=plot(WS,y2,'o','Color',blueShades(5,:),'MarkerSize',7,'LineWidth',1.5);
% plot3=plot(WS,y3,'o','Color',blueShades(4,:),'MarkerSize',7,'LineWidth',1.5);
% plot4=plot(WS,y4,'o','Color',blueShades(3,:),'MarkerSize',7,'LineWidth',1.5);
% plot5=plot(WS,y5,'o','Color',blueShades(2,:),'MarkerSize',7,'LineWidth',1.5);
% plot6=plot(WS,y6,'o','Color',blueShades(1,:),'MarkerSize',7,'LineWidth',1.5);

% % Perovskite panel
plot1=plot(WS,y1,'o','Color',redShades(6,:),'MarkerSize',7,'LineWidth',1.5);
hold on
plot2=plot(WS,y2,'o','Color',redShades(5,:),'MarkerSize',7,'LineWidth',1.5);
plot3=plot(WS,y3,'o','Color',redShades(4,:),'MarkerSize',7,'LineWidth',1.5);
plot4=plot(WS,y4,'o','Color',redShades(3,:),'MarkerSize',7,'LineWidth',1.5);
plot5=plot(WS,y5,'o','Color',redShades(2,:),'MarkerSize',7,'LineWidth',1.5);
plot6=plot(WS,y6,'o','Color',redShades(1,:),'MarkerSize',7,'LineWidth',1.5);

% lines
% plot7=plot(WS_0to4,y_0to4,':','Color',[0 0 0],'LineWidth',1.25);
% plot8=plot(WS_4to10,y_4to10,'--','Color',[0 0 0],'LineWidth',1.25);
plot9=plot(WS_0to10,y_0to10,'Color',[0 0 0],'LineWidth',1.25);

% % Sandia limits
% ylim([-5.1 -3.1])
% yticks(-5.1:0.2:-3.1)

% % Faiman limits
% ylim([20 110])
% yticks(20:10:110)

% Mattei, Skoplaki, and PVsyst limits
ylim([10 100])
yticks(10:10:100)

set(gca, 'FontSize', 17); 

hold off
set(gca, 'LineWidth', 1);
grid on;

xlabel('Wind speed (m/s)', FontSize=19);

% % Sandia y-label
% ylabel('log((T_m-T_a_m_b)/G)', FontSize=19);

% % Faiman y-label
% ylabel('G/(T_m-T_a_m_b)', FontSize=19);

% % Mattei y-label
% ylabel('^{G[(ατ)-η_r_e_f-βη_r_e_f(T_m-T_r_e_f)]}/_{(T_m-T_a_m_b)}', FontSize=19);

% Skoplaki y-label
% ylabel('^{(ατ)G[1-(η_{ref}/(ατ))]}/_{(T_c_e_l_l-T_a_m_b)}', FontSize=19);

% % PVsyst y-label
ylabel('^{G(ατ)(1-η_{ref})}/_{(T_{cell}-T_{amb})}', FontSize=19);
% % Sandia legend 1
% lgd1=legend([plot1 plot2 plot3 plot4 plot5 plot6],'200', '400', '600', '800', '1000', '1200', 'FontSize', 14, 'Location', 'northeast');

% Faiman, Mattei, Skoplaki, PVsyst legend 1
lgd1=legend([plot1 plot2 plot3 plot4 plot5 plot6],'200', '400', '600', '800', '1000', '1200', 'FontSize', 14, 'Location', 'northwest');
hold on;
title(lgd1,'Irr (W/m^2)')

ax1 = gca;
ax2 = axes('Position', get(ax1, 'Position'), 'Color', 'none');

% % Sandia legend 2
% lgd2=legend(ax2, [plot7 plot8 plot9], ['  0-4 m/s: y=',num2str(round(p_0to4(1,1),4)),'x',num2str(round(p_0to4(1,2),3))],['4-10 m/s: y=',num2str(round(p_4to10(1,1),4)),'x',num2str(round(p_4to10(1,2),3))], ['0-10 m/s: y=',num2str(round(p_0to10(1,1),4)),'x',num2str(round(p_0to10(1,2),3))], 'FontSize', 14, 'Location', 'southwest');

% Faiman, Mattei, Skoplaki, PVsyst legend 2
lgd2=legend(ax2, [plot9], ['0-10 m/s: y=',num2str(round(p_0to10(1,1),4)),'x+',num2str(round(p_0to10(1,2),3))], 'FontSize', 14, 'Location', 'southeast');

set(lgd2, 'Color', [1 1 1]);
set(ax2, 'Visible', 'off');
set(ax2, 'XTick', []);
set(ax2, 'YTick', []);

hold off;
