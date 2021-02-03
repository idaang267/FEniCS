clear all 
close all

Y = 2;
Nx = [250, 500, 625, 250, 500, 625, 250, 500, 625, 250, 500, 625];
Ny = [100, 200, 100, 100, 200, 100, 100, 200, 100, 100, 200, 100];
kappa = [10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10];
ellx = [5, 5, 5, 10, 10, 10, 15, 15, 15, 20, 20, 20];
colors = {[0, 0, 0], [0 0.4470 0.7410], [0.3010 0.7450 0.9330],...
          [0, 0, 0], [0 0.4470 0.7410], [0.3010 0.7450 0.9330],...
          [0, 0, 0], [0 0.4470 0.7410], [0.3010 0.7450 0.9330],...
          [0, 0, 0], [0 0.4470 0.7410], [0.3010 0.7450 0.9330]};

% colors = {[0, 0, 0], [0 0.4470 0.7410], [0.3010 0.7450 0.9330],[0.4940 0.1840 0.5560],...
%     [0, 0, 0], [0 0.4470 0.7410], [0.3010 0.7450 0.9330],[0.4940 0.1840 0.5560],...
%     [0, 0, 0], [0 0.4470 0.7410], [0.3010 0.7450 0.9330],[0.4940 0.1840 0.5560],...
%     [0, 0, 0], [0 0.4470 0.7410], [0.3010 0.7450 0.9330],[0.4940 0.1840 0.5560]};

SubPlotTitle = {"ellx = 5", "ellx = 10", "ellx = 15", "ellx = 20"};

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
    
    if i == 3 || i == 6 || i == 9
        count = count + 1; 
    end 
end

figure(2)
xlabel("Y-Coordinates")
ylabel("$P_{33}$", 'Interpreter', 'Latex')
leg = legend([h2(1), h2(2)], ["10","10000"], 'Location', 'Best');
title(leg, "$\kappa$", 'Interpreter', 'Latex');
saveas(gcf, "StressProfilekappa10.fig")

%%
xline(0, '--k')
for ii = 1:50
    xline(ii*hsize(end), '--k')
end
saveas(gcf, 'ResidualStressProfile.eps', 'epsc')