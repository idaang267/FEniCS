clear all
close all
y1_d = PureShearTestS1.y1;
y2_d = PureShearTestS1.y2;
y1_p_1 = PureShearTestS1.y3;
y2_p_1 = PureShearTestS1.y4;
y1_p_2 = PureShearTestS1.y5;
y2_p_2 = PureShearTestS1.y6;



%% Long and Hui displacement BCs
x = linspace(-3,0.0,500);
y = linspace(0.0,0.0,500);

r = sqrt(x.^2 + y.^2);
% Theta (rad) 
th = atan(y./x)+pi;
th(end) = pi; 
% h0 is the height 
h0 = 1.0;
% Amount of displacement 
Delta = 0.363366; 
% Stretch
lambda_a = 1 + Delta/h0;

% Known a for specific geometry 
a = 2*sqrt(h0/pi)*(lambda_a - 1/lambda_a);
% Unknown amplitudes
c = 1.25; %1.36;
% Displacement BC: 
y2 = a*sqrt(r).*sin(th./2); 
y1 = c*r.*cos(th);

figure(5)
hold on 
h1 = plot(y1,y2, 'ko', 'LineWidth', 2.5);
hold on 
h3 = plot(y1_p_1, y2_p_1, 'bo', 'LineWidth', 2.5);
h4 = plot(y1_p_2, y2_p_2, 'go', 'LineWidth', 2.5);
% h5 = plot(y1_p_3, y2_p_3, 'mo', 'LineWidth', 2.5);
% h6 = plot(y1_p_4, y2_p_4, 'ko', 'LineWidth', 2.5);

xlabel("$y_1$", 'Interpreter', 'LaTeX')
ylabel("$y_2$", 'Interpreter', 'LaTeX')
ax.FontSize = 15; 
leg = legend([h1, h3, h4],["Analytical", "X=5, Y=5", "pf=X=5, Y=2"], 'Location', 'Best');
% leg = legend([h2, h3, h4],["Discrete", "Original", "Threshold: 0.85"], 'Location', 'Best');
% leg = legend([h3, h4, h5, h6],["Original", "Threshold: 0.9", "Threshold: 0.8", "Threshold: 0.7"], 'Location', 'Best');
% leg = legend([h1, h3, h4, h5],["Analytical", "X=10, Y=2", "X=10, Y=5", "X=10, Y=10"], 'Location', 'Best');
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')

%%
y = linspace(0.0,0.0,500);

x_l = linspace(-3,0,500);
x_u = linspace(-3,2,500); 

r_l = sqrt(x_l.^2 + y.^2);
r_u = sqrt(x_u.^2 + y.^2);

% Theta (rad) 
th_l = atan(y./x_l)+pi;
th_l(end) = pi; 
th_u = atan(y./x_u)+pi;
th_u(end) = pi; 
% h0 is the height 
h0 = 1.0;
% Amount of displacement 
Delta = 0.363366; 
% Stretch
lambda_a = 1 + Delta/h0;

% Known a for specific geometry 
a = 2*sqrt(h0/pi)*(lambda_a - 1/lambda_a);
% Unknown amplitudes
c_low = 1.15;
c_upp = 1.555;

% Displacement BC: 
y2_l = a*sqrt(r_l).*sin(th_l./2); 
y1_l = c_low*r_l.*cos(th_l);

y2_u = a*sqrt(r_u).*sin(th_u./2);
y1_u = c_upp*r_u.*cos(th_u);

figure(2) 
hold on
a = area(y1_l,y2_l, 'FaceColor',[0 0.4470, 0.7410],'EdgeColor', 'none');
a.FaceAlpha = 0.5;
a1 = area(y1_u,y2_u, 'FaceColor', 'white', 'EdgeColor', 'none');
h1 = plot(y1_d,y2_d,'Marker','o', 'MarkerFaceColor',[0.4940 0.1840 0.5560],'Color',[0.4940 0.1840 0.5560]);
h2 = plot(y1_p_1,y2_p_1, 'Marker','o', 'MarkerFaceColor','k', 'Color', 'k');
leg = legend([a, h1, h2],["Asymptotic Solution", "Phase Field", "Discrete"], 'Location', 'Best');

xlabel("$y_1$", 'Interpreter', 'LaTeX')
ylabel("$y_2$", 'Interpreter', 'LaTeX')
ax.FontSize = 15; 
set(gca,'fontsize',15)
set(gca, 'Layer', 'top')

set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
xlim([-10,-0.001]);
ylim([0.01,1]);
saveas(gcf, 'Fig2B_trim.eps', 'epsc')
saveas(gcf, "Fig2B_trim.fig")
% xlim([-1,-0.01]);
% ylim([0.05,1]);
%%
y = linspace(0.0,0.0,500);

x = linspace(-8,3,500);
r = sqrt(x.^2 + y.^2);

% Theta (rad) 
th = atan(y./x)+pi;
th(end) = pi; 
% h0 is the height 
h0 = 1.0;
% Amount of displacement 
Delta = 0.363366; 
% Stretch
lambda_a = 1 + Delta/h0;

% Known a for specific geometry 
a = 2*sqrt(h0/pi)*(lambda_a - 1/lambda_a);
% Unknown amplitudes
c = 1.25;

% Displacement BC: 
y2 = a*sqrt(r).*sin(th./2);
y1 = c*r.*cos(th);

figure(2) 
hold on
a = plot(y1,y2, 'Color',[0 0.4470, 0.7410], 'LineWidth', 2);
h1 = plot(y1_d,y2_d,'Marker','o', 'MarkerFaceColor',[0.4940 0.1840 0.5560],'Color',[0.4940 0.1840 0.5560]);
h2 = plot(y1_p_1,y2_p_1, 'Marker','o', 'MarkerFaceColor','k', 'Color', 'k');
leg = legend([a, h1, h2],["Asymptotic Solution", "Phase Field", "Discrete"], 'Location', 'Best');

xlabel("$y_1$", 'Interpreter', 'LaTeX')
ylabel("$y_2$", 'Interpreter', 'LaTeX')
ax.FontSize = 15; 
set(gca,'fontsize',15)
set(gca, 'Layer', 'top')

set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
xlim([-10,-0.001]);
ylim([0.005,1]);
saveas(gcf, 'Fig2B_c1_25.eps', 'epsc')
saveas(gcf, "Fig2B_c1_25.fig")

%%
y = linspace(0.0,0.0,500);
x = linspace(-8,3,500);
r = sqrt(x.^2 + y.^2);

% Theta (rad) 
th = atan(y./x)+pi;
th(end) = pi; 
% h0 is the height 
h0 = 1.0;
% Amount of displacement 
Delta = 0.363366; 
% Stretch
lambda_a = 1 + Delta/h0;

% Known a for specific geometry 
a = 2*sqrt(h0/pi)*(lambda_a - 1/lambda_a);
% Unknown amplitudes
c = 1.25;

% Displacement BC: 
y2 = a*sqrt(r).*sin(th./2);
y1 = c*r.*cos(th);

figure(2) 
hold on
a = plot(y1,y2, 'Color',[0 0.4470, 0.7410], 'LineWidth', 2);
h1 = plot(y1_d(1:993),y2_d(1:993),'Marker','o', 'MarkerFaceColor',[0.4940 0.1840 0.5560],'Color',[0.4940 0.1840 0.5560]);
h2 = plot(y1_p_1(1:993),y2_p_1(1:993), 'Marker','o', 'MarkerFaceColor','k', 'Color', 'k');
leg = legend([a, h1, h2],["Asymptotic Solution", "Phase Field", "Discrete"], 'Location', 'Best');

xlabel("$y_1$", 'Interpreter', 'LaTeX')
ylabel("$y_2$", 'Interpreter', 'LaTeX')
ax.FontSize = 15; 
set(gca,'fontsize',15)
set(gca, 'Layer', 'top')

set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
xlim([-10,-0.01]);
ylim([0.05,1]);
saveas(gcf, 'Fig2B_trim.eps', 'epsc')
saveas(gcf, "Fig2B_trim.fig")