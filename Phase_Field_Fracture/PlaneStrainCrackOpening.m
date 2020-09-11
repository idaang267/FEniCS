clear all 
close all

% Figure 3 replication 

% Theta (rad) 
th = linspace(-pi, pi, 200);

% Choose x= b0/a^(2-2n) for Figure 4 

x = [0.5, 1, 5, 10];

n = [1, 1, 1, 1];

for ii = 1:length(x)
    % Kappa 
    k = 1 - 1/n(ii); 
    
    % omega [dimensionless function]
    w = sqrt(1 - (k*sin(th)).^2);

    % U (dimensionless) 
    U = sin(th./2).*sqrt(1 - (2*(k*cos(th./2)).^2)/(1+w)).*(w+k.*cos(th)).^(k/2);

    y2 = real(U); 
    y1 = -x(ii)*y2.^2; 
    
    figure(1)
    plot(y1, y2, 'LineWidth', 2.5)
    hold on
end 

%% Plotting

ax = gca;
xlabel("$y_1/a^{2n}$", 'Interpreter', 'LaTeX')
ylabel("$y_2/a^{2n}$", 'Interpreter', 'LaTeX')
xlim([-0.2 0])
ylim([0 0.4])
leg = legend(["0.5", "1", "5", "10"], 'Location', 'Best');
title(leg, "$b_0/a^{2-2n}$", 'Interpreter', 'Latex');
ax.FontSize = 20; 

saveas(gcf, 'Fig4A.eps', 'epsc')
% saveas(gcf, "Fig4.fig")