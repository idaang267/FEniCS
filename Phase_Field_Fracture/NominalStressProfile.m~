clear all 
close all

Y = 2;
Nx = [250, 250, 250, 250, 500, 500, 500, 500, 625, 625, 625, 625, 1250, 1250, 1250, 1250];
Ny = [100, 100, 100, 100, 200, 200, 200, 200, 100, 100, 100, 100, 100, 100, 100, 100];
kappa = [10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10];
ellx = [5, 10, 15, 20, 5, 10, 15, 20, 5, 10, 15, 20, 5, 10, 15, 20];
% colors = {[0, 0, 0], [0 0.4470 0.7410],...
%           [0, 0, 0], [0 0.4470 0.7410],...
%           [0, 0, 0], [0 0.4470 0.7410],...
%           [0, 0, 0], [0 0.4470 0.7410]};

colors = {[0, 0, 0], [0 0.4470 0.7410], [0.3010 0.7450 0.9330],[0.4940 0.1840 0.5560],...
    [0, 0, 0], [0 0.4470 0.7410], [0.3010 0.7450 0.9330],[0.4940 0.1840 0.5560],...
    [0, 0, 0], [0 0.4470 0.7410], [0.3010 0.7450 0.9330],[0.4940 0.1840 0.5560],...
    [0, 0, 0], [0 0.4470 0.7410], [0.3010 0.7450 0.9330],[0.4940 0.1840 0.5560]};

SubPlotTitle = {"hsize = 0.02", "hsize = 0.01", "hsize = 0.008", "hsize = 0.004"};
% SubPlotTitle = {"ellx = 5", "ellx = 10", "ellx = 15", "ellx = 20"};
% SubPlotTitle = {"Width: 0.40", "Width: 0.80", "Width = 0.16", "Width = 0.32"};

% DataPath = '/home/ida/Documents/GitHub/FEniCS/Phase_Field_Fracture/StressProfile/';
DataPath = '/Users/idaang/Documents/GitHub/FEniCS/Phase_Field_Fracture/StressProfile/';
count = 1; 
for i = 1:length(ellx)
    hsize(i) = Y/Ny(i);
    PhaseFieldWidth(i) = 4*hsize(i)*ellx(i);
    DataMatrix = importfile([DataPath 'Nx_' num2str(Nx(i)) '_Ny_' num2str(Ny(i)) '_ellx_' num2str(ellx(i)) '_k_' num2str(kappa(i)) '.csv' ] );
    
    Points(i,:) = DataMatrix.Points0; 
    NominalStress(i,:) = DataMatrix.NominalStress4; 
    
    figure(2)
    subplot(2,2,count)     
    title(SubPlotTitle{count})
    hold on 
    h2(i) = plot(Points(i,:), NominalStress(i,:), 'Color', colors{1,i}, 'LineWidth', 2.5);
    set(gca, 'YScale', 'log')
    ylim([1E-4 1])
    
    if i == 4 || i == 8 || i == 12
        count = count + 1; 
    end 
end

figure(2)
xlabel("Y-Coordinates")
ylabel("$P_{33}$", 'Interpreter', 'Latex')
leg = legend([h2(1), h2(2), h2(3), h2(4)], ["20: 0.4", "40: 0.8", "60: 1.2", "80: 1.6"], 'Location', 'Best');
leg = legend([h2(5), h2(6), h2(7), h2(8)], ["20: 0.2", "40: 0.4", "60: 0.6", "80: 0.8"], 'Location', 'Best');
leg = legend([h2(9), h2(10), h2(11), h2(12)], ["20: 0.16", "40: 0.32", "60: 0.48", "80: 0.64"], 'Location', 'Best');
leg = legend([h2(13), h2(14), h2(15), h2(16)], ["20: 0.08", "40: 0.16", "60: 0.24", "80: 0.32"], 'Location', 'Best');
title(leg, "Elements");
title(leg, "$\kappa$", 'Interpreter', 'Latex');
saveas(gcf, "StressProfile_Nx_v_Ny_v_k_10_hsize_v.fig")

%%
xline(0, '--k')
for ii = 1:50
    xline(ii*hsize(end), '--k')
end
saveas(gcf, 'ResidualStressProfile.eps', 'epsc')