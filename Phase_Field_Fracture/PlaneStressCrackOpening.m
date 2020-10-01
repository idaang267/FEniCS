clear all 
close all

% Theta (rad) 
th = linspace(-pi, pi, 200);
% h0 is the height 
h0 = 1.5;
% Amount of displacement 
Delta = 1.0; 
% stretch in direction of displacement 
lambda_a = 1 + Delta/h0;

% Known a for specific geometry 
ac = 2*sqrt(h0/pi)*(lambda_a - inv(lambda_a));
a = [ac];
% Unknown amplitudes
c = [1.55 1.15];

for ii = 1:length(c)
    y2(ii,:) = (1/a)*sin(th./2);
    y1(ii,:) = 2*c(ii)*y2(ii,:).^2*(cos(th)/(1-cos(th)));    
end 

%% Plotting
y1_1 = Compiled.y1_norm;
y2_1 = Compiled.y2_norm;
y1_2 = Compiled.y1_norm1;
y2_2 = Compiled.y2_norm1;
y1_3 = Compiled.y1_norm2;
y2_3 = Compiled.y2_norm2;

figure(1)
hold on 
plot(y1(1,:), y2(1,:), 'LineWidth', 2.5)
plot(y1(2,:), y2(2,:), 'LineWidth', 2.5)

h1 = plot(y1_1, y2_1, 'LineWidth', 2.5);
h2 = plot(y1_2, y2_2, 'LineWidth', 2.5);
h3 = plot(y1_3, y2_3, 'LineWidth', 2.5);
leg = legend([h1, h2, h3],["0.75+0.5l", "0.75+0.25l", "0.75+0.1l"], 'Location', 'Best');

ax = gca;
xlabel("$y_1/a^{2n}$", 'Interpreter', 'LaTeX')
ylabel("$y_2/a^{2n}$", 'Interpreter', 'LaTeX')
xlim([-0.4 0])
ylim([0 0.4])
ax.FontSize = 20; 
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')

saveas(gcf, 'PlaneStressCrackOpening.eps', 'epsc')
% saveas(gcf, "Fig4.fig")

% title(leg, "$c$", 'Interpreter', 'Latex');