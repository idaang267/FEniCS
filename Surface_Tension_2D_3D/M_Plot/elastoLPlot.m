clear all 
close all

% directory = ["Case_chi_0.1_g_0.00_l0_2.50_exp_1.1/",...
%              "Case_chi_0.1_g_0.10_l0_2.50_exp_1.1/",...
%              "Case_chi_0.1_g_0.20_l0_2.50_exp_1.1/",...
%              "Case_chi_0.1_g_0.30_l0_2.50_exp_1.1/",...
%              "Case_chi_0.1_g_0.40_l0_2.50_exp_1.1/",...
%              "Case_chi_0.1_g_0.50_l0_2.50_exp_1.1/",...
%              "Case_chi_0.1_g_0.60_l0_2.50_exp_1.1/",...
%              "Case_chi_0.1_g_0.70_l0_2.50_exp_1.1/",...
%              "Case_chi_0.1_g_0.80_l0_2.50_exp_1.1/",...
%              "Case_chi_0.1_g_0.90_l0_2.50_exp_1.1/",...
%              "Case_chi_0.1_g_1.00_l0_2.50_exp_1.1/"];
% 
% chi = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1];
% n = 10^-3;
% l0 = [2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5];
% J0 = l0.^3; 
% gamma = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0];

% directory = ["exp_chi_0.2_g_0.00_l0_2.50_exp_1.1/",...
%              "exp_chi_0.2_g_0.10_l0_2.50_exp_1.1/",...
%              "exp_chi_0.2_g_0.20_l0_2.50_exp_1.1/",...
%              "exp_chi_0.2_g_0.30_l0_2.50_exp_1.1/",...
%              "exp_chi_0.2_g_0.40_l0_2.50_exp_1.1/",...
%              "exp_chi_0.2_g_0.50_l0_2.50_exp_1.1/",...
%              "exp_chi_0.2_g_0.60_l0_2.50_exp_1.1/",...
%              "exp_chi_0.2_g_0.70_l0_2.50_exp_1.1/",...
%              "exp_chi_0.2_g_0.80_l0_2.50_exp_1.1/",...
%              "exp_chi_0.2_g_0.90_l0_2.50_exp_1.1/",...
%              "exp_chi_0.2_g_1.00_l0_2.50_exp_1.1/"];
%      
% chi = [0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2];
% n = 10^-3;
% l0 = [2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5];
% J0 = l0.^3; 
% gamma = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0];

directory = ["Case_chi_0.3_g_0.00_l0_2.50_exp_1.1/",...
             "Case_chi_0.3_g_0.10_l0_2.50_exp_1.1/",...
             "Case_chi_0.3_g_0.20_l0_2.50_exp_1.1/",...
             "Case_chi_0.3_g_0.30_l0_2.50_exp_1.1/",...
             "Case_chi_0.3_g_0.40_l0_2.50_exp_1.1/",...
             "Case_chi_0.3_g_0.50_l0_2.50_exp_1.1/",...
             "Case_chi_0.3_g_0.60_l0_2.50_exp_1.1/",...
             "Case_chi_0.3_g_0.70_l0_2.50_exp_1.1/",...
             "Case_chi_0.3_g_0.80_l0_2.50_exp_1.1/",...
             "Case_chi_0.3_g_0.90_l0_2.50_exp_1.1/",...
             "Case_chi_0.3_g_1.00_l0_2.50_exp_1.1/"];

chi = [0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3];
n = 10^-3;
l0 = [2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5];
J0 = l0.^3; 
gamma = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0];

% directory = ["Case_chi_0.4_g_0.00_l0_2.50_exp_1.1/",...
%              "Case_chi_0.4_g_0.10_l0_2.50_exp_1.1/",...
%              "Case_chi_0.4_g_0.20_l0_2.50_exp_1.1/",...
%              "Case_chi_0.4_g_0.30_l0_2.50_exp_1.1/",...
%              "Case_chi_0.4_g_0.40_l0_2.50_exp_1.1/",...
%              "Case_chi_0.4_g_0.50_l0_2.50_exp_1.1/",...
%              "Case_chi_0.4_g_0.60_l0_2.50_exp_1.1/",...
%              "Case_chi_0.4_g_0.70_l0_2.50_exp_1.1/",...
%              "Case_chi_0.4_g_0.80_l0_2.50_exp_1.1/",...
%              "Case_chi_0.4_g_0.90_l0_2.50_exp_1.1/",...
%              "Case_chi_0.4_g_1.00_l0_2.50_exp_1.1/"];
%          
% % These calculations happen automatically for legend generation. 
% chi = [0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4];
% n = 10^-3;
% l0 = [2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5];
% J0 = l0.^3; 
% gamma = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0];
% 
% directory = ["Case_chi_0.5_g_0.00_l0_2.50_exp_1.1/",...
%              "Case_chi_0.5_g_0.10_l0_2.50_exp_1.1/",...
%              "Case_chi_0.5_g_0.20_l0_2.50_exp_1.1/",...
%              "Case_chi_0.5_g_0.30_l0_2.50_exp_1.1/",...
%              "Case_chi_0.5_g_0.40_l0_2.50_exp_1.1/",...
%              "Case_chi_0.5_g_0.50_l0_2.50_exp_1.1/",...
%              "Case_chi_0.5_g_0.60_l0_2.50_exp_1.1/",...
%              "Case_chi_0.5_g_0.70_l0_2.50_exp_1.1/",...
%              "Case_chi_0.5_g_0.80_l0_2.50_exp_1.1/",...
%              "Case_chi_0.5_g_0.90_l0_2.50_exp_1.1/",...
%              "Case_chi_0.5_g_1.00_l0_2.50_exp_1.1/"];
%          
% % These calculations happen automatically for legend generation. 
% chi = [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5];
% n = 10^-3;
% l0 = [2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5];
% J0 = l0.^3; 
% gamma = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0];

% Elastocapillary Lengthscale 
le = gamma.*l0;

for i = 1:length(directory)
    % Import data including time, number of steps, and radius
    dataplot = importfile0(directory(i) + "data_plot.txt", 1, inf);
    time = dataplot.VarName1;
    steps = dataplot.VarName2;
    % Non-Normalized radius 
    radius = dataplot.VarName3;
    % Normalized radius by initial radius 
    radNorm = radius./radius(1);
        
    minMaxRad(i) = radius(1);
    minMaxRad(i+length(directory)) = radius(end); 
    
    minMaxRadNorm(i) = radNorm(1);
    minMaxRadNorm(i+length(directory)) = radNorm(end); 
    
    %Calculate the point at which the swelling radius is 99% to equilibrium
    % value. Must calculate start point after ramping 
    % Non-Normalized 
    radDiff(i) = radius(end) - radius(6); 
    rad99(i) = radius(end) - radDiff(i)*0.01;
    % Normalized 
    radDiffNorm(i) = radNorm(end) - radNorm(6);
    radNorm99(i) = radNorm(end) - radDiffNorm(i)*0.01; 
    
    % Compute the closest value in a vector, radius, to 99% of swelling 
    % radius. Find both value and index in array radius 
    [minValue(i),closeIndex(i)] = min(abs(radius-rad99(i)));
    % Find time value at the particular index
    timeValue(i) = time(closeIndex(i));
    % Find the radius value at the particular index 
    radValue(i) = radius(closeIndex(i)); 
    radNormValue(i) = radNorm(closeIndex(i));
    
    %% Figure 
    minRadVal = min(minMaxRad);
    maxRadVal = max(minMaxRad);
    
    figure(2)
    h1(i) = plot(le(i), timeValue(i), '*', 'LineWidth', 1.5, 'MarkerSize', 15);
    hold on
    %ylim([1000, 10000])
    xlabel("l", 'Interpreter', 'LaTeX')
    ylabel("$t^{99}_{eqm}$", 'Interpreter', 'LaTeX')
    set(gca,'fontsize',22)
    
    leg{i} = num2str(gamma(i), '%.2f\n');
    
end

%%
fo2 = fitoptions('Method','NonlinearLeastSquares',...
               'Robust','Off',...
               'Algorithm', 'Trust-Region',...
               'MaxIter', 400);
           
[curve2,gof2] = fit(le', timeValue','exp2', fo2)

%%
a1 = 4785;
b1 = -1.655;
c1 = 3227;
d1 = -0.2595;
x = linspace(0,2.5,1000); 
y1 = a1*exp(b1*x) + c1*exp(d1*x);

a2 = 2907;
b2 = -2.18;
c2 = 4388;
d2 = -0.4546;
y2 = a2*exp(b2*x) + c2*exp(d2*x);

a3 = 2662; 
b3 = -2.114;
c3 = 3366;
d3 = -0.386;
y3 = a3*exp(b3*x) + c3*exp(d3*x);

a4 = 1220;
b4 = -3.654;
c4 = 3325;
d4 = -0.4935;
y4 = a4*exp(b4*x) + c4*exp(d4*x);

a5 = 885.5;
b5 = -1.955;
c5 = 1460;
d5 = -0.3037;
y5 = a5*exp(b5*x) + c5*exp(d5*x);

figure(3)
hold on 
h1 = plot(x, y1, 'LineWidth', 2.5);
h2 = plot(x, y2, 'LineWidth', 2.5);
h3 = plot(x, y3, 'LineWidth', 2.5);
%h0 = scatter(le, timeValue, 'fill');
hold off
xlabel("$\hat{l}$", 'Interpreter', 'LaTeX')
ylabel("$\hat{t}^{99}_{eqm}$", 'Interpreter', 'LaTeX')
set(gca,'fontsize',22)
leg = legend([h1 h2 h3], {"0.1", "0.2", "0.3"}, 'Location', 'Best');
title(leg, {"$\chi$"}, 'Interpreter', 'Latex');

saveas(gcf, 'bp_chi_04_l0_2_radius_norm.eps', 'epsc')
saveas(gcf, "bp_chi_04_l0_2_radius_norm.fig")


