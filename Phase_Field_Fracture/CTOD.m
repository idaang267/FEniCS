clear all 
close all

% Theta (rad) 
th = linspace(-pi, pi, 501);

% Choose n = 1 for neo-hookean model

% Fig 3A 
x = [0.5 1 5 10];
% Kappa 
k = 1 - 1/1;
% omega [dimensionless function]
w = sqrt(1 - (k*sin(th)).^2);

for ii = 1:length(x)
    
    % U (dimensionless) 
    U = sin(th./2).*sqrt(1 - (2*(k*cos(th./2)).^2)/(1+w)).*(w+k.*cos(th)).^(k/2);

    y2 = real(U);
    y1 = -x(ii)*y2.^2;
    figure(1)
    plot(y1,y2, 'LineWidth', 2.5)
    hold on 
end 

xlabel("$y_1 / a^{2n}$", 'Interpreter', 'LaTeX')
ylabel("$y_2 / a^{2n}$", 'Interpreter', 'LaTeX')
ax = gca; 
ax.FontSize = 15; 
xlim([-0.2 0])
ylim([0 0.4])
leg = legend(["0.5", "1", "5", "10"], 'Location', 'Best');
title(leg, "$b_0/a^{2-2n}$", 'Interpreter', 'Latex');
saveas(gcf, 'Fig4.eps', 'epsc')
