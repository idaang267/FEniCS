%% Plotting for Figure 2 in Paper

%% Read in file 
Disp_PF = DEnergiesS2.Displacement; 
EE_PF = DEnergiesS2.ElasticEnergy; 
DE_PF = DEnergiesS2.DissipatedEnergy;
Tot_PF = DEnergiesS2.TotalEnergy;

Disp_Dis = DEnergiesS2.Displacement1; 
EE_Dis = DEnergiesS2.ElasticEnergy1; 
DE_Dis = DEnergiesS2.DissipatedEnergy1;
Tot_Dis = DEnergiesS2.TotalEnergy1;

% x = 0.4659*ones(1,100);
x = 0.5454*ones(1,100); 
y = linspace(0,7,100);

x_v = linspace(0, 0.55, 100);
y_v = 3*ones(1,100);

%%
figure
hold on 
h1 = plot(Disp_PF,EE_PF,'Color', [0 0.4470 0.7410],  'LineWidth', 2.5)
h2 = plot(Disp_PF,DE_PF, 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 2.5)
h3 = plot(Disp_PF,Tot_PF, 'Color', [0.4940 0.1840 0.5560], 'LineWidth', 3)
plot(x,y, 'k:', 'LineWidth', 2.5)
% plot(x_v,y_v, 'k:', 'LineWidth', 2.5)

xlabel("Displacement")
ylabel("Energy")
leg = legend([h1, h2, h3],"Elastic", "Dissipated", "Total")
set(gca,'fontsize',20)
xlim([0 .55])

saveas(gcf, 'Fig2A_gce.eps', 'epsc')
saveas(gcf, "Fig2A_gce.fig")
%%
figure
hold on 
h4 = plot(Disp_Dis,EE_Dis, 'Color', [0 0.4470 0.7410], 'LineWidth', 3)
h5 = plot(Disp_Dis,DE_Dis, 'Color', [0.8500 0.3250 0.0980],'LineWidth', 3)
h6 = plot(Disp_Dis,Tot_Dis, 'Color', [0.4940 0.1840 0.5560], 'LineWidth', 3)
plot(x,y, 'k:', 'LineWidth', 2.5)
% plot(x_v,y_v, 'k:', 'LineWidth', 2.5)

xlabel("Displacement")
ylabel("Energy")
leg = legend([h4, h5, h6],"Elastic", "Dissipated", "Total")
set(gca,'fontsize',20)
xlim([0 .55])

saveas(gcf, 'Fig2B_gce.eps', 'epsc')
saveas(gcf, "Fig2B_gce.fig")
