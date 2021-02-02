clear all 
close all

Y = 2;
Nx = [250, 250, 250, 250, 250, 250, 250, 250,...
     250, 250, 250, 250, 250, 250, 250, 250];
Ny = [100, 100, 100, 100, 100, 100, 100, 100,...
     100, 100, 100, 100, 100, 100, 100, 100];
kappa = [10, 100, 1000, 10000, 10,100,1000,10000,...
         10, 100, 1000, 10000, 10, 100, 1000, 10000];
ellx = [5, 5, 5, 5, 10, 10, 10, 10,...
        15, 15, 15, 15, 20, 20, 20, 20];
    
colors = {[0, 0, 0], [0 0.4470 0.7410], [0.3010 0.7450 0.9330],[0.4940 0.1840 0.5560],...
    [0, 0, 0], [0 0.4470 0.7410], [0.3010 0.7450 0.9330],[0.4940 0.1840 0.5560],...
    [0, 0, 0], [0 0.4470 0.7410], [0.3010 0.7450 0.9330],[0.4940 0.1840 0.5560],...
    [0, 0, 0], [0 0.4470 0.7410], [0.3010 0.7450 0.9330],[0.4940 0.1840 0.5560]};

SubPlotTitle = {"ellx = 5", "ellx = 10", "ellx = 15", "ellx = 20"};

% DataPath = '/home/ida/Documents/GitHub/FEniCS/Phase_Field_Fracture/Test_ellx_change/';
DataPath = '/Users/idaang/Documents/GitHub/FEniCS/Phase_Field_Fracture/Test_ellx_change/';
count = 1; 
for i = 1:length(ellx)
    hsize(i) = Y/Ny(i);
    PhaseFieldWidth(i) = 4*hsize(i)*ellx(i);
    DataMatrix = importfile([DataPath 'Nx_' num2str(Nx(i)) '_Ny_' num2str(Ny(i)) '_ellx_' num2str(ellx(i)) '_k_' num2str(kappa(i)) '.csv' ] );
    
    Points(i,:) = DataMatrix.Points1; 
    Pressure(i,:) = DataMatrix.Pressure;
    NominalStress(i,:) = DataMatrix.NominalStress4; 
    
    figure(1)
    subplot(2,2,count)
    title(SubPlotTitle{count})
    if i == 5 || i == 9 || i == 13 
        count = count + 1; 
    end 
    hold on 
    h1(i) = plot(Points(i,:), Pressure(i,:), 'Color', colors{1,i}, 'LineWidth', 2.5);
    xlim([-0.2 0.2])

    figure(2)
    subplot(2,2,count)
    title(SubPlotTitle{count})
    hold on 
    h2(i) = plot(Points(i,:), NominalStress(i,:), 'Color', colors{1,i}, 'LineWidth', 2.5);

end

figure(1)
xlabel("Y-Coordinates")
ylabel("Pressure")

xline(0, '--k')
for ii = 1:20
    xline(ii*hsize(end), '--k')
end
% leg = legend(["5: " + num2str(PhaseFieldWidth(1)), "5: " + num2str(PhaseFieldWidth(2)), "5: " + num2str(PhaseFieldWidth(3)), "5: " + + num2str(PhaseFieldWidth(4))], 'Location', 'Best');
leg = legend(["10", "100", "1000", "10000"], 'Location', 'Best');
title(leg, "$\kappa$", 'Interpreter', 'Latex');
saveas(gcf, "PressureProfile.fig")

figure(2) 
xlabel("Y-Coordinates")
ylabel("$P_{33}$", 'Interpreter', 'Latex')
saveas(gcf, "StressProfile.fig")

xline(0, '--k')
for ii = 1:50
    xline(ii*hsize(end), '--k')
end
% leg = legend(["5: " + num2str(PhaseFieldWidth(1)), "5: " + num2str(PhaseFieldWidth(2)), "5: " + num2str(PhaseFieldWidth(3)), "5: " + + num2str(PhaseFieldWidth(4))], 'Location', 'Best');
leg = legend(["10", "100", "1000", "10000"], 'Location', 'Best');
title(leg, "$\kappa$", 'Interpreter', 'Latex');
saveas(gcf, 'ResidualStressProfile.eps', 'epsc')
