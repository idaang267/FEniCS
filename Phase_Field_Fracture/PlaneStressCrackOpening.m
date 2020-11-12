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
    y2(ii,:) = a*sin(th./2);
    y1(ii,:) = -c(ii)*(y2(ii,:)/a).^2;
    
    y2_norm(ii,:) = (1/a)*sin(th./2);
    y1_norm(ii,:) = 2*c(ii)*y2(ii,:).^2*(cos(th)/(1-cos(th)));    
end 

%% Plotting
% y1_1 = Compiled.y1_norm + .005;
y1_1 = TrimS1.y1;
y2_1 = TrimS1.y2;
y1_2 = TrimS1.y3;
y2_2 = TrimS1.y4;

y1_1 = Case1S1.y1;
y2_1 = Case1S1.y2;
y1_2 = Case1S1.y3;
y2_2 = Case1S1.y4;
y1_3 = Case1S1.y5;
y2_3 = Case1S1.y6;
y1_4 = Case1S1.y7;
y2_4 = Case1S1.y8;

y1_1 = Case2S1.y1;
y2_1 = Case2S1.y2;
y1_2 = Case2S1.y3;
y2_2 = Case2S1.y4;
y1_3 = Case2S1.y5;
y2_3 = Case2S1.y6;
y1_4 = Case2S1.y7;
y2_4 = Case2S1.y8;

y1_1 = Case3S1.y1;
y2_1 = Case3S1.y2;
y1_2 = Case3S1.y3;
y2_2 = Case3S1.y4;
y1_3 = Case3S1.y5;
y2_3 = Case3S1.y6;
y1_4 = Case3S1.y7;
y2_4 = Case3S1.y8;

figure(1)
hold on 
plot(-y1(1,:), y2(1,:), 'LineWidth', 2.5)
plot(-y1(2,:), y2(2,:), 'LineWidth', 2.5)
hold on
h1 = plot(y1_1, y2_1, 'o');
h2 = plot(y1_2, y2_2, 'LineWidth', 2.5);
h3 = plot(y1_3, y2_3, 'LineWidth', 2.5);
h4 = plot(y1_4, y2_4, 'LineWidth', 2.5);
h5 = plot(y1_5, y2_5, 'LineWidth', 2.5);
h6 = plot(y1_6, y2_6, 'LineWidth', 2.5);

title(leg, "$No Phase-Field$", 'Interpreter', 'Latex');
leg = legend([h1, h2],["X+2, Y+5", "Trimmed"], 'Location', 'Best');
leg = legend([h1, h2, h3, h4],["X+4, Y+5", "X+3, Y+6", "X+2, Y+8", "X+1, Y+10"], 'Location', 'Best');
leg = legend([h1, h2, h3, h4],["X+2, Y+5", "X+4, Y+5", "X+6, Y+5", "X+8, Y+5"], 'Location', 'Best');
leg = legend([h1, h2, h3, h4],["X+3, Y+5", "X+3, Y+6", "X+3, Y+8", "X+3, Y+10"], 'Location', 'Best');

% leg = legend([h1, h2, h3, h4],["Discrete Stabilized", "Phase-Field Stabilized", "Displacement", "Taylor-Hood"], 'Location', 'Best');

ax = gca;
xlabel("$y_1$", 'Interpreter', 'LaTeX')
ylabel("$y_2$", 'Interpreter', 'LaTeX')
xlim([0 0.5])
ylim([0 1.0])
ax.FontSize = 15; 
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')

saveas(gcf, 'Send.eps', 'epsc')

%%
figure(2)
Points = Vertical.Points1;
a = Vertical.a;
Damage = Vertical.Damage;
h0 = plot(Points, Damage, 'LineWidth', 2.5)
hold on
h1 = plot(Points, a, 'LineWidth', 2.5)
leg = legend([h0, h1],["Damage", "a"], 'Location', 'Best');
