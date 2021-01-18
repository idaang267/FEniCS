clear all 
close all
ux_f = TestS2.uxux00;
uy_f = TestS2.uyuy00;

x_f = TestS2.Points0;
y1_f = TestS2.y1;
y2_f = TestS2.y2; 

%% Long and Hui displacement BCs
x = linspace(-2.5,0,200);
y = linspace(0.0,0.0,200);

r = sqrt(x.^2 + y.^2);
% Theta (rad) 
th = atan(y./x)+pi;
th(end) = pi; 
% h0 is the height 
h0 = 1.0;
% Amount of displacement 
Delta = 0.01; 
% Stretch
lambda_a = 1 + Delta/h0;

% Known a for specific geometry 
a = 2*sqrt(h0/pi)*(lambda_a - 1/lambda_a);
% Unknown amplitudes
c = 1.55 ;
% Displacement BC: 
y2 = a*sqrt(r).*sin(th./2); 
y1 = c*r.*cos(th);

figure(5)
h1 = plot(y1,y2, 'ro', 'LineWidth', 2.5)
hold on 
h2 = plot(y1_f, y2_f, 'bo-', 'LineWidth', 2.5)
xlabel("$y_1$", 'Interpreter', 'LaTeX')
ylabel("$y_2$", 'Interpreter', 'LaTeX')
ax.FontSize = 15; 
leg = legend([h1, h2],["Analytical", "FEA"], 'Location', 'Best');

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

xlim([0 0.04])
ylim([0 0.04])
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
saveas(gcf, 'Log.eps', 'epsc')

%     y1(ii,:) = 2*c(ii)*((y2(ii,:)/a).^2).*(cos(th)/(1-cos(th)));
%     y1(ii,:) = c(ii)*(y2(ii,:)/a).^2;
    
%     y2_norm(ii,:) = (1/a)*sin(th./2);
%     y1_norm(ii,:) = 2*c(ii)*y2(ii,:).^2*(cos(th)/(1-cos(th)));    

%% Submodel 
r = 0.1;
x_coord = 0.08;
y_coord = sqrt(r^2 - x_coord^2)

x = linspace(x_coord,y_coord,200);
y = linspace(0.0,0.0,200);

r = sqrt(x.^2 + y.^2);
% Theta (rad) 
th = atan(y./x)+pi;
th(end) = pi; 
% h0 is the height 
h0 = 1.0;
% Amount of displacement 
Delta = 0.01; 
% Stretch
lambda_a = 1 + Delta/h0;
% Known a for specific geometry 
a = 2*sqrt(h0/pi)*(lambda_a - 1/lambda_a);
% Unknown amplitudes
c = 1.55 ;
% Displacement BC: 
y2 = a*sqrt(r).*sin(th./2); 
y1 = c*r.*cos(th);

figure(5)
hold on 
plot(x,y1, 'r', 'LineWidth', 2.5)
plot(x,y2, 'b', 'LineWidth', 2.5)

y1_f = TestS1.y1;
y2_f = TestS1.y2; 

figure(6)
hold on
plot(y1,y2, 'k', 'LineWidth', 2.5) 
plot(y1_f, y2_f, 'b', 'LineWidth', 2.5)
title(leg, "$c$", 'Interpreter', 'Latex');
leg = legend([h1,h2,h3,h4,h5,h6,h7,h8,h9],["1.15", "1.20", "1.25", "1.30", "1.35", "1.40", "1.45", "1.50", "1.55"], 'Location', 'Best');


title(leg, "$No Phase-Field$", 'Interpreter', 'Latex');
leg = legend([h1, h2],["X+2, Y+5", "Trimmed"], 'Location', 'Best');
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

saveas(gcf, 'Varying_c.eps', 'epsc')

%% 

% th = linspace(-pi, pi, 200);

u1 = 1.15*r.*cos(-th);
u2 = 1.4741*sqrt(r).*sin(-th./2);
figure(1)
hold on 
plot(x, u1, 'LineWidth', 2.5)
plot(x, u2, 'LineWidth', 2.5)

u1_test = Untitled.Displacement0;
u2_test = Untitled.Displacement1;
x_test = Untitled.Points0;



figure(2)
hold on
plot(x_test, u1_test, 'LineWidth', 2.5);
plot(x_test, u2_test, 'LineWidth', 2.5);
