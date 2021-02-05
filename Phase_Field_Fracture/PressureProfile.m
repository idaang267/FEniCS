clear all 
close all

Y = 2;
Nx = [1250, 1250, 1250, 1250, 1250, 1250, 1250, 1250];
Ny = [100, 100, 100, 100, 100, 100, 100, 100];
kappa = [10, 10, 10, 10, 10000, 10000, 10000, 10000];
ellx = [5, 10, 15, 20, 5, 10, 15, 20];
% colors = {[0, 0, 0], [0 0.4470 0.7410],...
%           [0, 0, 0], [0 0.4470 0.7410],...
%           [0, 0, 0], [0 0.4470 0.7410],...
%           [0, 0, 0], [0 0.4470 0.7410]};

colors = {[0, 0, 0], [0 0.4470 0.7410], [0.3010 0.7450 0.9330],[0.4940 0.1840 0.5560],...
    [0, 0, 0], [0 0.4470 0.7410], [0.3010 0.7450 0.9330],[0.4940 0.1840 0.5560],...
    [0, 0, 0], [0 0.4470 0.7410], [0.3010 0.7450 0.9330],[0.4940 0.1840 0.5560],...
    [0, 0, 0], [0 0.4470 0.7410], [0.3010 0.7450 0.9330],[0.4940 0.1840 0.5560]};

SubPlotTitle = {"kappa: 10", "kappa: 10000"};
% SubPlotTitle = {"Width: 0.40", "Width: 0.8", "Width = 0.16", "Width = 0.32"};
% SubPlotTitle = {"ellx = 5", "ellx = 10", "ellx = 15", "ellx = 20"};

% DataPath = '/home/ida/Documents/GitHub/FEniCS/Phase_Field_Fracture/PressureProfile/';
DataPath = '/Users/idaang/Documents/GitHub/FEniCS/Phase_Field_Fracture/PressureProfile/';
count = 1; 
for i = 1:length(ellx)
    hsize(i) = Y/Ny(i);
    PhaseFieldWidth(i) = 4*hsize(i)*ellx(i);
    DataMatrix = importfile([DataPath 'Nx_' num2str(Nx(i)) '_Ny_' num2str(Ny(i)) '_ellx_' num2str(ellx(i)) '_k_' num2str(kappa(i)) '.csv' ] );
    
    Points(i,:) = DataMatrix.Points1; 
    Pressure(i,:) = DataMatrix.Pressure;
    
    figure(1)
        subplot(1,2,count)
%     subplot(2,2,count)
    title(SubPlotTitle{count})
    if i == 4 || i == 8% || i == 6
        count = count + 1; 
    end 
    hold on 
    h1(i) = plot(Points(i,:), Pressure(i,:), 'Color', colors{1,i}, 'LineWidth', 2.5);
    xlim([-0.25 0.25])
%     ylim([-8 0])
end

figure(1)
xlabel("Y-Coordinates")
ylabel("Pressure")
xlim([-0.05, +0.05])

% leg = legend(["5: " + num2str(PhaseFieldWidth(1)), "5: " + num2str(PhaseFieldWidth(2)), "5: " + num2str(PhaseFieldWidth(3)), "5: " + + num2str(PhaseFieldWidth(4))], 'Location', 'Best');
leg = legend([h1(1), h1(2)], ["10", "10000"], 'Location', 'Best');
leg = legend([h1(1), h1(2), h1(3), h1(4)], ["20", "40", "60", "80"], 'Location', 'Best');
leg = legend([h1(1), h1(2)], ["0.02", "0.01"], 'Location', 'Best');
leg = legend([h1(3), h1(4)], ["0.02", "0.01"], 'Location', 'Best');
leg = legend([h1(5), h1(6)], ["0.008", "0.004"], 'Location', 'Best');
leg = legend([h1(7), h1(8)], ["0.008", "0.004"], 'Location', 'Best');
title(leg, "Elements");

title(leg, "$\kappa$", 'Interpreter', 'Latex');
saveas(gcf, "PressureProfile_Nx_1250_Ny_100.fig")
