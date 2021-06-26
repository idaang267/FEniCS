
%% 

Displacement = double(stabilizedenergies.e00);
ElasticEnergy = double(stabilizedenergies.e15);
SurfaceEnergy = double(stabilizedenergies.e1);
TotalEnergy = ElasticEnergy + SurfaceEnergy; 

hold on 
h1 = plot(Displacement(1:end-2),ElasticEnergy(1:end-2), 'LineWidth', 2.5)
h2 = plot(Displacement(1:end-2),SurfaceEnergy(1:end-2), 'LineWidth', 2.5)
h3 = plot(Displacement(1:end-2),TotalEnergy(1:end-2), 'LineWidth', 2.5)
leg = legend([h1, h2, h3],"Elastic", "Dissipated", "Total")
xlabel("Displacement")
ylabel("Energies")

xlim([0,Displacement(end-2)])
set(gca,'fontsize',15)
xlim([0,3]);

saveas(gcf, 'Fig3C.eps', 'epsc')
saveas(gcf, "Fig3C.fig")

%% 
Disp_PF = DEnergiesS2.Displacement; 
EE_PF = DEnergiesS2.ElasticEnergy; 
DE_PF = DEnergiesS2.DissipatedEnergy;
Tot_PF = DEnergiesS2.TotalEnergy;

Disp_Dis = DEnergiesS2.Displacement1; 
EE_Dis = DEnergiesS2.ElasticEnergy1; 
DE_Dis = DEnergiesS2.DissipatedEnergy1;
Tot_Dis = DEnergiesS2.TotalEnergy1;

y = linspace(0,7,100);
x = 0.5454*ones(1,100);

hold on 
h1 = plot(Disp_PF,EE_PF,'Color', [0 0.4470 0.7410],  'LineWidth', 3)
h2 = plot(Disp_PF,DE_PF, 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 3)
h3 = plot(Disp_PF,Tot_PF, 'Color', [0.4940 0.1840 0.5560], 'LineWidth', 3)

h4 = plot(Disp_Dis,EE_Dis, 'Color', [0 0.4470 0.7410],'LineStyle', '-.', 'LineWidth', 3.5)
h5 = plot(Disp_Dis,DE_Dis, 'Color', [0.8500 0.3250 0.0980], 'LineStyle', '-.','LineWidth', 3.5)
h6 = plot(Disp_Dis,Tot_Dis, 'Color', [0.4940 0.1840 0.5560], 'LineStyle', '-.', 'LineWidth', 3.5)

plot(x,y, 'k:', 'LineWidth', 2.5)

xlabel("Displacement")
ylabel("Energies")
leg = legend([h1, h2, h3],"Elastic", "Dissipated", "Total")
set(gca,'fontsize',20)

%%
figure
hold on 
h1 = plot(Disp_PF,EE_PF,'Color', [0 0.4470 0.7410],  'LineWidth', 2.5)
h2 = plot(Disp_PF,DE_PF, 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 2.5)
h3 = plot(Disp_PF,Tot_PF, 'Color', [0.4940 0.1840 0.5560], 'LineWidth', 3)
plot(x,y, 'k:', 'LineWidth', 2.5)

xlabel("Displacement")
ylabel("Energy")
leg = legend([h1, h2, h3],"Elastic", "Dissipated", "Total")
set(gca,'fontsize',20)
xlim([0 .6])

%%
figure
hold on 
h4 = plot(Disp_Dis,EE_Dis, 'Color', [0 0.4470 0.7410], 'LineWidth', 3)
h5 = plot(Disp_Dis,DE_Dis, 'Color', [0.8500 0.3250 0.0980],'LineWidth', 3)
h6 = plot(Disp_Dis,Tot_Dis, 'Color', [0.4940 0.1840 0.5560], 'LineWidth', 3)
plot(x,y, 'k:', 'LineWidth', 2.5)

xlabel("Displacement")
ylabel("Energy")
leg = legend([h4, h5, h6],"Elastic", "Dissipated", "Total")
set(gca,'fontsize',20)
xlim([0 .6])

