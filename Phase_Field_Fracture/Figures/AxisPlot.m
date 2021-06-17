Points = DJDataS2.Points0;
J002 = DJDataS2.Y002; 
J012 = DJDataS2.Y012;
D002 = DJDataS2.Y1;
D012 = DJDataS2.Y112; 

%%
hold on 
yyaxis left 
h1 = plot(Points, J002, 'LineWidth', 2.5)
h2 = plot(Points, J012, 'LineWidth', 2.5)
% plot(Points, ones(length(Points)), 'k:')
ylabel('J')
yyaxis right
h3 = plot(Points, D002, 'LineWidth', 2.5) 
h4 = plot(Points, D012, 'LineWidth', 2.5) 
legend([h1, h3, h2, h4], "Y = 0.02", "Y = 0.02", "Y = 0.12", "Y = 0.12")
ylim([-0.01, 1])
ylabel('Damage')
xlabel('X-Coordinates')
ax = gca; 
ax.FontSize = 20;
xlim([0, 0.25])

%%
Disp = DEnergiesS2.Displacement; 
PFEla = DEnergiesS2.ElasticEnergy; 
PFDis = DEnergiesS2.DissipatedEnergy;
PFTot = DEnergiesS2.TotalEnergy;

DisEla = DEnergiesS2.ElasticEnergy1; 
DisDis = DEnergiesS2.DissipatedEnergy1;
DisTot = DEnergiesS2.TotalEnergy1;

Analytical = 0.4866*ones(10);
Y = linspace(0,7,10); 

hold on 
h1 = plot(Disp, PFEla, 'Color', '#0072BD', 'LineWidth', 2.5)
h2 = plot(Disp, PFDis,'Color', '#D95319', 'LineWidth', 2.5)
h3 = plot(Disp, PFTot,'Color', '#7E2F8E', 'LineWidth', 2.5)

h4 = plot(Disp, DisEla,'Color', '#0072BD', 'LineStyle', ":", 'LineWidth', 2.5)
h5 = plot(Disp, DisDis, 'Color', '#D95319', 'LineStyle', ":", 'LineWidth', 2.5)
h6 = plot(Disp, DisTot, 'Color', '#7E2F8E', 'LineStyle', ":", 'LineWidth', 2.5)

h7 = plot(Analytical, Y, 'k:', 'LineWidth', 2)
xlabel('Displacement')
ylabel('Energies')
legend([h1, h2, h3], "Elastic",  "Dissipated", "Total")
legend([h4, h5, h6], "Elastic",  "Dissipated", "Total")
ax = gca; 
ax.FontSize = 20;
xlim([0, 0.55])
