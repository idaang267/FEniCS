clear all 
close all

y1_d_ls = PhaseFieldPlotS2.u1;
y2_d_ls = PhaseFieldPlotS2.u2; 
y1_p_1 = PhaseFieldPlotS2.u3;
y2_p_1 = PhaseFieldPlotS2.u4; 
y1_p_3 = TestS1.y7;
y2_p_3 = TestS1.y8; 

y1_d_ls = PhaseFieldPlotS1.y1;
y2_d_ls = PhaseFieldPlotS1.y2; 
y1_p_1 = PhaseFieldPlotS1.y3;
y2_p_1 = PhaseFieldPlotS1.y4;
y1_p_2 = PhaseFieldPlotS1.y5;
y2_p_2 = PhaseFieldPlotS1.y6;
y1_p_3 = PhaseFieldPlotS1.y7;
y2_p_3 = PhaseFieldPlotS1.y8;

y1_d = PhaseFieldS1.y1;
y2_d = PhaseFieldS1.y2;
y1_p_1 = PhaseFieldS1.y3;
y2_p_1 = PhaseFieldS1.y4;

y1_d = PhaseFieldS1.y5;
y2_d = PhaseFieldS1.y6;
y1_p_1 = PhaseFieldS1.y7;
y2_p_1 = PhaseFieldS1.y8;

%% Long and Hui displacement BCs
x = linspace(-4,0.00,500);
y = linspace(0.00,0.00,500);

r = sqrt(x.^2 + y.^2);
% Theta (rad) 
th = atan(y./x)+pi;
th(end) = pi; 
% h0 is the height 
h0 = 1.0;
% Amount of displacement 
Delta = 0.1; 
% Stretch
lambda_a = 1 + Delta/h0;

% Known a for specific geometry 
a = 2*sqrt(h0/pi)*(lambda_a - 1/lambda_a);
% Unknown amplitudes
c = 2.1; %1.36;
% Displacement BC: 
y2 = a*sqrt(r).*sin(th./2); 
y1 = c*r.*cos(th);

figure(5)
hold on 
h1 = plot(y1,y2, 'ko', 'LineWidth', 2.5);
h2 = plot(y1_d, y2_d, 'ro', 'LineWidth', 2.5)
hold on 
h3 = plot(y1_p_1, y2_p_1, 'bo', 'LineWidth', 2.5);
h4 = plot(y1_p_2, y2_p_2, 'go', 'LineWidth', 2.5);
h5 = plot(y1_p_3, y2_p_3, 'mo', 'LineWidth', 2.5);
xlabel("$y_1$", 'Interpreter', 'LaTeX')
ylabel("$y_2$", 'Interpreter', 'LaTeX')
ax.FontSize = 15; 
leg = legend([h2, h3],["Discrete", "Phase Field"], 'Location', 'Best');
% leg = legend([h1, h3, h4, h5],["Analytical", "X=10, Y=2", "X=10, Y=5", "X=10, Y=10"], 'Location', 'Best');
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')

figure(6) 
hold on 
plot(x, y2)
plot(x_p, y2_p_1)


%%
saveas(gcf, "LargeStrain_h_0-001.fig")
title("$\lambda_a = 1.1$", "Interpreter", 'LaTeX')
leg = legend([h1, h3, h4],["Analytical", "Phase Field: kappa = 10E3", "Phase Field: kappa = 10E4"], 'Location', 'Best');

figure(6)
h1 = plot(y1,y2, 'ro', 'LineWidth', 2.5)
hold on 
h2 = plot(y1_d_ss, y2_d_ss, 'ko', 'LineWidth', 2.5)
h3 = plot(y1_p_ss, y2_p_ss, 'bo', 'LineWidth', 2.5)
xlabel("$y_1$", 'Interpreter', 'LaTeX')
ylabel("$y_2$", 'Interpreter', 'LaTeX')
ax.FontSize = 15; 
leg = legend([h1, h2, h3],["Analytical", "Discrete", "Phase-Field"], 'Location', 'Best');
xlim([-0.2 0])
ylim([-2 0])
set(gca, 'XScale', 'linear')
set(gca, 'YScale', 'linear')
saveas(gcf, "PhaseFieldComparison_3.fig")
saveas(gcf, 'Vary_D0.01.eps', 'epsc')

figure(1)
h1 = plot(x, y1(1,:), 'r', 'LineWidth', 2.5)
hold on 
h2 = plot(x_f, y1_f, 'b', 'LineWidth', 2.5)

figure(2)
hold on
plot(x, y2(1,:), 'r', 'LineWidth', 2.5)
plot(x_f, y2_f, 'b', 'LineWidth', 2.5)

h1 = plot(-y1(1,:), y2(1,:), 'r', 'LineWidth', 2.5);
hold on 
h2 = plot(y1_f, y2_f, 'b', 'LineWidth', 2.5);

% h1 = plot(y1_1, y2_1, 'o');



%     y1(ii,:) = 2*c(ii)*((y2(ii,:)/a).^2).*(cos(th)/(1-cos(th)));
%     y1(ii,:) = c(ii)*(y2(ii,:)/a).^2;
    
%     y2_norm(ii,:) = (1/a)*sin(th./2);
%     y1_norm(ii,:) = 2*c(ii)*y2(ii,:).^2*(cos(th)/(1-cos(th)));    

%% Submodel 
x = linspace(-0.5,0,200);
y = linspace(0.0,0.0,200);
r = sqrt(x.^2 + y.^2);
th = atan(y./x) + pi;% Theta (rad) 
th(end) = pi; 
plot(x,th) 

% h0 is the height 
h0 = 1.0;
% Amount of displacement 
Delta = 0.01; 
% Stretch
lambda_a = 1 + Delta/h0;
% Known a for specific geometry 
a = 2*sqrt(h0/pi)*(lambda_a - 1/lambda_a);
% Unknown amplitudes
c = 2.15;
% Displacement BC: 
y2 = a*sqrt(r).*sin(th./2); 
y1 = c*r.*cos(th);

% Import 
u1_f = SubModelS1.Displacement0;
u2_f = SubModelS1.Displacement1; 
x_f = SubModelS1.Points0;
y_f = SubModelS1.Points1;
y1_f = SubModelS1.y1;
y2_f = SubModelS1.y2; 

figure(3)
hold on 
plot(x,y1, 'k', 'LineWidth', 2.5)
plot(x_f,u1_f, 'r--', 'LineWidth', 2.5)
leg = legend("Asymptotic","FEA Displacement", 'Location', 'Best');
xlabel("X-Coordinates", 'Interpreter', 'LaTeX')
ylabel("X-Displacement", 'Interpreter', 'LaTeX')

figure(4)
hold on 
plot(x,y2, 'k', 'LineWidth', 2.5)
plot(x_f,u2_f, 'r', 'LineWidth', 2.5)
leg = legend("Asymptotic","FEA Displacement", 'Location', 'Best');
xlabel("X-Coordinates", 'Interpreter', 'LaTeX')
ylabel("Y-Displacement", 'Interpreter', 'LaTeX')

figure(6)
hold on
plot(y1,y2, 'k', 'LineWidth', 2.5) 
plot(u1, u2, 'b', 'LineWidth', 2.5)


title(leg, "$c$", 'Interpreter', 'Latex');
leg = legend([h1,h2,h3,h4,h5,h6,h7,h8,h9],["1.15", "1.20", "1.25", "1.30", "1.35", "1.40", "1.45", "1.50", "1.55"], 'Location', 'Best');


title(leg, "$No Phase-Field$", 'Interpreter', 'Latex');
leg = legend([h1, h2, h3, h4],["X+4, Y+5", "X+3, Y+6", "X+2, Y+8", "X+1, Y+10"], 'Location', 'Best');
leg = legend([h1, h2, h3, h4],["X+2, Y+5", "X+4, Y+5", "X+6, Y+5", "X+8, Y+5"], 'Location', 'Best');
leg = legend([h1, h2, h3, h4],["X+3, Y+5", "X+3, Y+6", "X+3, Y+8", "X+3, Y+10"], 'Location', 'Best');

% leg = legend([h1, h2, h3, h4],["Discrete Stabilized", "Phase-Field Stabilized", "Displacement", "Taylor-Hood"], 'Location', 'Best');

ax = gca;
xlabel("$y_1$", 'Interpreter', 'LaTeX')
ylabel("$y_2$", 'Interpreter', 'LaTeX')
xlim([0 0.2])
ylim([0 0.2])
ax.FontSize = 15; 
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')

% saveas(gcf, 'Varying_c.eps', 'epsc')

%%
h = 1;
mu = 1;
Gc = 6.033135859170121141e+06;

syms lambda_a
eq = Gc == h*mu*(lambda_a - 1/lambda_a)^2;
ans = double(solve(eq, lambda_a));
