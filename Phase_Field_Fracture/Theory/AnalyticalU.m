clear all 
close all

% Figure 3 replication 

% Theta (rad) 
th = linspace(-pi, pi, 101);

% Choose n = 1 for neo-hookean model

% Fig 3A 
n = [0.7 1 2 3];
% Fig 3B
%n = [1.7, 2, 3];

for ii = 1:length(n)
    % Kappa 
    k = 1 - 1/n(ii); 
    
    % omega [dimensionless function]
    w = sqrt(1 - (k*sin(th)).^2);

    % U (dimensionless) 
    U = sin(th./2).*sqrt(1 - (2*(k*cos(th./2)).^2)/(1+w)).*(w+k.*cos(th)).^(k/2);

    figure(1) 
    plot(th, U, 'LineWidth', 2.5)
    hold on 
    
end 

%% Plotting
figure(1) 
grid on 
xlabel("$\theta$", 'Interpreter', 'LaTeX')
ylabel("$U(\theta, n)$", 'Interpreter', 'LaTeX')
ax = gca; 
ax.XTick = -pi:pi/2:pi;
ax.XTickLabel = {'-\pi','-\pi/2','0','\pi/2','\pi'};
ax.FontSize = 15; 
xlim([-pi pi])
leg = legend(["n = 0.7", "n = 1", "n = 2", "n = 3"], 'Location', 'Best');
saveas(gcf, 'Fig3A.eps', 'epsc')