

Displacement = double(stabilizedenergies.e00);
ElasticEnergy = double(stabilizedenergies.e16);
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

saveas(gcf, 'Fig2C.eps', 'epsc')
saveas(gcf, "Fig2C.fig")