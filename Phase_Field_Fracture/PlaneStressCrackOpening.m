clear all 
close all
x = linspace(-0.5,0,200);
y = linspace(0.005,0.005,200);

r = sqrt(x.^2 + y.^2);

% Theta (rad) 
th = linspace(0, pi, 200);
% h0 is the height 
h0 = 10000.0;
% Amount of displacement 
Delta = 100; 
% stretch in direction of displacement 
lambda_a = 1 + Delta/h0;

% Known a for specific geometry 
ac = 2*sqrt(h0/pi)*(lambda_a - inv(lambda_a));
a = [ac];
% Unknown amplitudes
c = [1.55 1.15];
% c = [0.422449];
for ii = 1:length(c)
    y2(ii,:) = a*sin(th./2);
    y1(ii,:) = 2*c(ii)*((y2(ii,:)/a).^2).*(cos(th)/(1-cos(th)));
%     y1(ii,:) = c(ii)*cos(th);
%     y1(ii,:) = c(ii)*(y2(ii,:)/a).^2;
    
    y2_norm(ii,:) = (1/a)*sin(th./2);
    y1_norm(ii,:) = 2*c(ii)*y2(ii,:).^2*(cos(th)/(1-cos(th)));    
end

figure(1)
hold on 
plot(-y1(1,:), y2(1,:), 'LineWidth', 2.5)
plot(-y1(2,:), y2(2,:), 'LineWidth', 2.5)
hold on
% h1 = plot(y1_1, y2_1, 'o');

y1_1 = Disp.y5;
y2_1 = Disp.y6; 

h1 = plot(y1_1, y2_1, 'LineWidth', 2.5);
xlabel("$y_1$", 'Interpreter', 'LaTeX')
ylabel("$y_2$", 'Interpreter', 'LaTeX')
ax.FontSize = 15; 
xlim([0 0.04])
ylim([0 0.04])
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
saveas(gcf, 'Log.eps', 'epsc')

%% Plotting

y1_1 = VaryingC.y1;
y2_1 = VaryingC.y2;
y1_2 = VaryingC.y3;
y2_2 = VaryingC.y4;
y1_3 = VaryingC.y5;
y2_3 = VaryingC.y6;
y1_4 = VaryingC.y7;
y2_4 = VaryingC.y8;
y1_5 = VaryingC.y9;
y2_5 = VaryingC.y10;
y1_6 = VaryingC.y11;
y2_6 = VaryingC.y12;
y1_7 = VaryingC.y13;
y2_7 = VaryingC.y14;
y1_8 = VaryingC.y15;
y2_8 = VaryingC.y16;
y1_9 = VaryingC.y17;
y2_9 = VaryingC.y18;

figure(1)
hold on
h1 = plot(y1_1, y2_1, 'LineWidth', 2.5);
h2 = plot(y1_2, y2_2, 'LineWidth', 2.5);
h3 = plot(y1_3, y2_3, 'LineWidth', 2.5);
h4 = plot(y1_4, y2_4, 'LineWidth', 2.5);
h5 = plot(y1_5, y2_5, 'LineWidth', 2.5);
h6 = plot(y1_6, y2_6, 'LineWidth', 2.5);
h7 = plot(y1_7, y2_7, 'LineWidth', 2.5);
h8 = plot(y1_8, y2_8, 'LineWidth', 2.5);
h9 = plot(y1_9, y2_9, 'LineWidth', 2.5);
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
xlim([0 0.5])
ylim([0 0.2])
ax.FontSize = 15; 
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')

saveas(gcf, 'Varying_c_comparison.eps', 'epsc')

%%
figure(2)
Points = Vertical.Points1;
a = Vertical.a;
Damage = Vertical.Damage;
h0 = plot(Points, Damage, 'LineWidth', 2.5)
hold on
h1 = plot(Points, a, 'LineWidth', 2.5)
leg = legend([h0, h1],["Damage", "a"], 'Location', 'Best');

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
