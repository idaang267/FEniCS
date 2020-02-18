clear all
close all

% Comparison of chi = 0.2 with and without surface tension
% directory = ["exp_chi_0.2_g_1.00_l0_2.00_exp_1.6/", "exp_chi_0.2_g_1.00_l0_2.50_exp_1.6/", "exp_chi_0.2_g_1.00_l0_3.00_exp_1.6/",...
%              "exp_chi_0.2_g_0.00_l0_2.00_exp_1.6/", "exp_chi_0.2_g_0.00_l0_2.50_exp_1.6/", "exp_chi_0.2_g_0.00_l0_3.00_exp_1.6/"];

% Comparison of chi = 0.4 with and without surface tension
% directory = ["exp_chi_0.4_g_1.00_l0_2.00_exp_1.6/", "exp_chi_0.4_g_1.00_l0_2.50_exp_1.6/", "exp_chi_0.4_g_1.00_l0_3.00_exp_1.6/",...
%              "exp_chi_0.4_g_0.00_l0_2.00_exp_1.6/", "exp_chi_0.4_g_0.00_l0_2.50_exp_1.6/", "exp_chi_0.4_g_0.00_l0_3.00_exp_1.6/"];

% Comparison of varying chi with surface tension
% directory = ["exp_chi_0.2_g_1.00_l0_2.00_exp_1.6/", "exp_chi_0.2_g_1.00_l0_2.50_exp_1.6/", "exp_chi_0.2_g_1.00_l0_3.00_exp_1.6/",...
%              "exp_chi_0.4_g_1.00_l0_2.00_exp_1.6/", "exp_chi_0.4_g_1.00_l0_2.50_exp_1.6/", "exp_chi_0.4_g_1.00_l0_3.00_exp_1.6/"];

directory = ["exp_chi_0.2_g_1.00_l0_2.00_exp_1.6/",...
             "alt_chi_0.2_lm_1.0_mu_1.0_l0_2.00/",...
             "alt_chi_0.2_lm_2.0_mu_2.0_l0_2.00/",...
             "alt_chi_0.2_lm_10_mu_10_l0_2.00/"];

         %% Parameters to control the figure legend 
% Controls the number of directories before change 
types = 1; 
dashType = {'-', '--', '--'};
% Controls asterik size
markSize = [5, 5, 5];
% Controls color of lines using rgb values
% colors = {[0.06 0.08 0.8], [0 0.5 0.9], [0.3 0.7 0.9],...
%           [0.8500 0.3250 0.0980], [1 0 0], [0.6350 0.0780 0.1840], [0 0 0]};

colors = {[1 0 0], [0, 0, 1], [0.06 0.08 0.8], [0.3 0.7 0.9], [0 0 0]};

% These calculations of chemical potential happen automatically for legend
% generation. 
n = 10^-3;      % N Omega
% With or without surface tension 
gamma = [1.0, 0.0, 0.0];
chi = [0.2, 0.2, 0.2];
l0 = [2.0, 2.0, 2.0];
J0 = l0.^3;
mu0 = n.*(1./l0-1./J0) + 1./J0 + log(1 - 1./J0) + chi./(J0.^2)+ n.*gamma*2./l0;

% Loop through Directories
for i = 1:length(directory)
    % Import data including time, number of steps, and radius
    dataplot = importfile(directory(i) + "data_plot.txt",1,inf);
    time = dataplot.e00;
    steps = dataplot.e1;
    % Non-Normalized radius 
    radius = dataplot.e2;
    % Normalized radius by initial radius 
    radNorm = radius./radius(1);
    
    minMaxRad(i) = radius(1);
    minMaxRad(i+length(directory)) = radius(end); 
    
    minMaxRadNorm(i) = radNorm(1);
    minMaxRadNorm(i+length(directory)) = radNorm(end); 
    
    % Calculate the point at which the swelling radius is 99% to equilibrium
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

    figure(1)
    grid on 
    set(gca,'fontsize',15, 'XScale', 'log')
    xlabel("Time")
    ylabel("Radius")
    h1(i) = plot(time,radius, dashType{i}, 'Color', colors{1,i}, 'LineWidth', 2);
    hold on
    h2(i) = plot(timeValue(i), radValue(i), '*', 'Color', colors{1,i}, 'MarkerSize', markSize(i));
    hold on 
    % Limits of x and y axis 
    xlim([1E-5 1E5]) 
    newLim = get(gca,'XLim'); 
    %newx = [1E-5, 1, 1E5];    
    newx = [1E-5, 1E-4, 1E-3, 1E-2,1E-1, 1, 1E1, 1E2, 1E3, 1E4,1E5];
    set(gca,'XTick', newx);     
    minRadVal = min(minMaxRad);
    maxRadVal = max(minMaxRad);
    ylim([minRadVal - 0.1, maxRadVal + 0.1])
    
    figure(2) 
    grid on
    set(gca,'fontsize',15, 'XScale', 'log')
    xlabel("Time")
    ylabel("R/R_{ini}")
    h3(i) = plot(time, radNorm, dashType{i}, 'Color', colors{1,i}, 'LineWidth', 2);
    hold on
    h4(i)= plot(timeValue(i), radNormValue(i), '*', 'Color', colors{1,i}, 'MarkerSize', markSize(i));
    hold on 
    % Limits of x and y axis 
    xlim([1E-5 1E5]) 
    newLim = get(gca,'XLim');
%    newx = [1E-5, 1, 1E5];        
    newx = [1E-5, 1E-4, 1E-3, 1E-2,1E-1, 1, 1E1, 1E2, 1E3, 1E4,1E5];
    set(gca,'XTick', newx); 

    minRadValNorm = min(minMaxRadNorm);
    maxRadValNorm = max(minMaxRadNorm);
    ylim([minRadValNorm - 0.1, maxRadValNorm + 0.1])
    
    % Legend 
    legStretch{i} = num2str(l0(i),'%.1f\n');
    legChem{i} = num2str(mu0(i),'%10.0e\n');
end 

figure(1) 
for ii = 1:types:length(directory)
    avgRad99 = mean(rad99(1,ii:ii+types-1));
    avgTimeVal = mean(timeValue(1,ii:ii+types-1));
    
    xV = avgTimeVal*ones(1,100); 
    yV = linspace(minRadVal - 0.1,maxRadVal + 0.1,100);
    plot(xV, yV, ':', 'Color', colors{1,end}, 'Linewidth', 2)
    hold on 
end 

% legend([h1(1), h1(2), h1(3), h1(4), h1(5), h1(6)],...
%         {"\lambda_0 = " + legStretch{1} +  ", \mu_0 = " + legChem{1},...
%          "\lambda_0 = " + legStretch{2} +  ", \mu_0 = " + legChem{2},...
%          "\lambda_0 = " + legStretch{3} +  ", \mu_0 = " + legChem{3},...
%          "\lambda_0 = " + legStretch{4} +  ", \mu_0 = " + legChem{4},...
%          "\lambda_0 = " + legStretch{5} +  ", \mu_0 = " + legChem{5},...
%          "\lambda_0 = " + legStretch{6} +  ", \mu_0 = " + legChem{6}
%          },...
%          'Location', 'Best')
legend([h1(1), h1(2), h1(3)],...
        {"\gamma = 1.0, \mu_0 = " + legChem{1},...
         "\lambda_s = \mu_s = 1.0, \mu_0 = " + legChem{2},...
         "\lambda_s = \mu_s = 2.0, \mu_0 = " + legChem{3}},...
         'Location', 'Best')     
saveas(gcf, "radius.pdf")
saveas(gcf, "radius.fig")

figure(2) 
for ii = 1:types:length(directory)
    avgRadNorm99 = mean(radNorm99(1,ii:ii+types-1));
    avgTimeVal = mean(timeValue(1,ii:ii+types-1));
    
    xV = avgTimeVal*ones(1,100); 
    yV = linspace(minRadValNorm - 0.1, maxRadValNorm + 0.1,100);
    plot(xV, yV, ':', 'Color', colors{1,end}, 'Linewidth', 2)
    hold on     
end 

% legend([h3(1), h3(2), h3(3), h3(4), h3(5), h3(6)],...
%         {"\lambda_0 = " + legStretch{1} +  ", \mu_0 = " + legChem{1},...
%          "\lambda_0 = " + legStretch{2} +  ", \mu_0 = " + legChem{2},...
%          "\lambda_0 = " + legStretch{3} +  ", \mu_0 = " + legChem{3},...
%          "\lambda_0 = " + legStretch{4} +  ", \mu_0 = " + legChem{4},...
%          "\lambda_0 = " + legStretch{5} +  ", \mu_0 = " + legChem{5},...
%          "\lambda_0 = " + legStretch{6} +  ", \mu_0 = " + legChem{6}
%          },...
%          'Location', 'Best')
legend([h3(1), h3(2), h3(3)],...
        {"\gamma = 1.0, \mu_0 = " + legChem{1},...
         "\lambda_s = \mu_s = 1.0, \mu_0 = " + legChem{2},...
         "\lambda_s = \mu_s = 2.0, \mu_0 = " + legChem{3}},...
         'Location', 'Best') 

saveas(gcf,'radius_norm.pdf')
saveas(gcf, "radius_norm.fig")


            
%%