clear all 
close all

% Theta (rad) 
th = linspace(-pi, pi, 200);
% h0 is the height 
h0 = 1.0;
% Amount of displacement 
Delta = 1.0; 
% stretch in direction of displacement 
lambda_a = 1 + Delta/h0;

% Known a for specific geometry 
ac = 2*sqrt(h0/pi)*(lambda_a - inv(lambda_a));
a = [ac];
% Unknown amplitudes
c = [2];

for ii = 1:length(a)
    y2 = (1/a(ii))*sin(th./2);
    y1 = 2*c(ii)*y2.^2*(cos(th)/(1-cos(th)));    
    figure(1)
    plot(y1, y2, 'LineWidth', 2.5)
    hold on
end 

%% Plotting

ax = gca;
xlabel("$y_1/a^{2n}$", 'Interpreter', 'LaTeX')
ylabel("$y_2/a^{2n}$", 'Interpreter', 'LaTeX')
xlim([-0.4 0])
ylim([0 0.4])
leg = legend(["0.5", "1", "1.15", "1.55"], 'Location', 'Best');
title(leg, "$c$", 'Interpreter', 'Latex');
ax.FontSize = 20; 

saveas(gcf, 'PlaneStressCrackOpening.eps', 'epsc')
% saveas(gcf, "Fig4.fig")