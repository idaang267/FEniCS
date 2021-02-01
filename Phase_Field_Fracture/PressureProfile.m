clear all 
close all

Y = 2;
Nx = 500;
Ny = 200;
ellx = [5,8,10,20];

hsize = Y/Ny;
limitX = ellx(end)*hsize;

DataPath = '/home/ida/Documents/GitHub/FEniCS/Phase_Field_Fracture/Test_ellx_change/';

for i = 1:length(ellx)
    PhaseFieldWidth(i) = 4*hsize*ellx(i);
    DataMatrix = importfile([DataPath 'Nx_' num2str(Nx) '_Ny_' num2str(Ny) '_ellx_' num2str(ellx(i)) '_k_10000.csv' ] );
    
    Points(i,:) = DataMatrix.Points1; 
    Pressure(i,:) = DataMatrix.Pressure;
    NominalStress(i,:) = DataMatrix.NominalStress4; 
    
    figure(1)
    hold on 
    plot(Points(i,:), Pressure(i,:), 'LineWidth', 2.5)
    
    figure(2)
    hold on 
    plot(Points(i,:), NominalStress(i,:), 'LineWidth', 2.5)
end

figure(1)
xlim([-limitX/2 limitX/2])
xlabel("Y-Coordinates")
ylabel("Pressure")

xline(0, '--k')
for ii = 1:(limitX/2)/hsize
    xline(ii*hsize, '--k')
end
leg = legend(["5: 0.2", "8: 0.32", "10: 0.4", "20: 0.8"], 'Location', 'Best');
title(leg, "$\ell$", 'Interpreter', 'Latex');
saveas(gcf, 'PressureProfile.eps', 'epsc')

figure(2) 
xlabel("Y-Coordinates")
ylabel("$P_{33}$", 'Interpreter', 'Latex')

xline(0, '--k')
for ii = 1:50
    xline(ii*hsize, '--k')
end
leg = legend(["5: 0.2", "8: 0.32", "10: 0.4", "20: 0.8"], 'Location', 'Best');
title(leg, "$\ell$", 'Interpreter', 'Latex');
saveas(gcf, 'ResidualStressProfile.eps', 'epsc')
