clear all
close all


% % Black, Purple, light blue, Black for dashed line 
% colors = {[0 0 0], [0 0 1], [0.3010 0.7450 0.9330],...
%           [0 0 0], [0 0 1], [0.3010 0.7450 0.9330],...
%           [0 0 0]};

% Black, Bright Blue, Med Blue, light blue, Purple, Black 
% colors = {[0, 0, 0],...
%         [0 0 1], [0 0.4470 0.7410], [0.3010 0.7450 0.9330],...
%         [0.4940 0.1840 0.5560], [0 0 0]};
    
%Red, Black, Bright Blue, Med Blue, Light Blue 
% colors = {[1 0 0], [0 0 0],...
%           [0 0 1], [0.3010 0.7450 0.9330],...
%           [0 0 0]};

% Black, Light blue, Black for dashed line 
% colors = {[0 0 0], [0 0.4470 0.7410],...
%           [0 0 0], [0 0.4470 0.7410], [0 0 0]};

% Three blues 
%colors = {[0 0 1], [0 0.4470 0.7410], [0.3010 0.7450 0.9330], [0 0 0]};

% Comparison of chi = 0.2 with and without surface tension
% directory = ["exp_chi_0.2_g_1.00_l0_2.00_exp_1.6/",...
%              "exp_chi_0.2_g_1.00_l0_2.50_exp_1.6/",...
%              "exp_chi_0.2_g_1.00_l0_3.00_exp_1.6/",...
%              "exp_chi_0.2_g_0.00_l0_2.00_exp_1.6/",...
%              "exp_chi_0.2_g_0.00_l0_2.50_exp_1.6/",...
%              "exp_chi_0.2_g_0.00_l0_3.00_exp_1.6/"];

% Comparison of varying chi with surface tension and fixed stretch
% directory = ["exp_chi_0.1_g_1.00_l0_2.00_exp_1.6/",...
%              "exp_chi_0.2_g_1.00_l0_2.00_exp_1.6/",...
%              "exp_chi_0.3_g_1.00_l0_2.00_exp_1.6/",...
%              "exp_chi_0.4_g_1.00_l0_2.00_exp_1.6/",...
%              "exp_chi_0.5_g_1.00_l0_2.00_exp_1.6/"];

% Comparison of Neo-Hookean vs Surface tension boundary potential 
% directory = ["exp_chi_0.4_g_0.00_l0_2.00_exp_1.6/",...
%              "fixed_chi_0.4_g_1.00_l0_2.00_exp_1.6/"...
%              "alt_fixed_chi_0.4_lm_10_mu_10_l0_2.0/",...
%              "alt_fixed_chi_0.4_lm_20_mu_20_l0_2.0/"];

% Comparison of chi = 0.2 with lower surface tension 
% directory = ["exp_chi_0.2_g_0.00_l0_2.50_exp_1.1/",...
%              "exp_chi_0.2_g_0.10_l0_2.50_exp_1.1/",...
%              "exp_chi_0.2_g_0.20_l0_2.50_exp_1.1/",...
%              "exp_chi_0.2_g_0.50_l0_2.50_exp_1.1/",...
%              "exp_chi_0.2_g_1.00_l0_2.50_exp_1.1/"];
      

% directory = ["C2_chi_0.4_g_0.00_l0_2.60/",...
%              "C1_chi_0.4_g_0.01_l0_2.60/",...
%              "C4_chi_0.4_g_0.00_l0_2.75/",...
%              "C3_chi_0.4_g_0.01_l0_2.75/"]; 
%          
% cases =  ["Pore_Case2_f_dict.mat",...
%           "Pore_Case1_f_dict.mat",...
%           "Pore_Case4_f_dict.mat",...
%           "Pore_Case3_f_dict.mat"];       

directory = ["exp_chi_0.2_l0_2/"]

%% Parameters to control the figure legend 
% Controls the number of directories before change 
types = 1; 
levPrec = 5; 
dashType = {'--', '-', '-.', '-.'};

% Controls asterik size
markSize = [15, 15, 15, 15, 15, 15, 15];

% Default MATLAB Colors, Blue, Orange, Purple, Light Blue, Dark Red, Green, 
colors = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.4940 0.1840 0.5560],...
         [0.3010 0.7450 0.9330], [0.6350 0.0780 0.1840], [0.4660 0.6740 0.1880],...
         [0 0 0]}; 
ConCol = 0; %length(cases); 

% These calculations of chemical potential happen automatically for legend
% generation. 
n = 10^-3;      % N Omega
% Second values are for linear results 
gamma = [1.00];
chi = [0.2];
l0 = [2];

% gamma = [0.00, 0.00, 0.0, 0.0];
% chi = [0.4, 0.4, 0.4, 0.4];
% l0 = [2, 2, 2, 2];
J0 = l0.^3;
mu0 = n.*(1./l0-1./J0) + 1./J0 + log(1 - 1./J0) + chi./(J0.^2)+ n.*gamma*2./l0;

for l = 1:length(gamma)
    % Legend 
    legChi{l} = num2str(chi(l), '%.1f\n');
    legStretch{l} = num2str(l0(l),'%.2f\n');
    legChem{l} = num2str(mu0(l),'%10.0e\n');
    legGamma{l} = num2str(gamma(l), '%.2f\n'); 
end 

% Automatically set lower/upper limits for the y-axis
for dir = 1:length(directory)+length(cases) 
    if dir <= length(directory) 
        % Import data just for radius 
        dataplot = importfile0(directory(dir) + "data_plot.txt", 1, inf);
        load(cases(dir))
    end
    % Non-Normalized radius
    radius = dataplot.VarName3;      
    % Normalize radius by initial radius 
    radNorm = radius./radius(1);        
    % Automatically setting upper and lower bounds for:
    % Non-normalized
    MinRad(dir) = radius(1);
    MaxRad(dir) = radius(end); 
    % Normalized 
    MinRadNorm(dir) = radNorm(1);
    MaxRadNorm(dir) = radNorm(end); 
    
    % Take into account the limits of the linear cases 
    if dir > length(directory) 
        % Reference and Current 
        fRef = (f+1)*l0(end);
        fCurr = f+1; 
        % Reference
        MinRad(dir) = fRef(1);
        MaxRad(dir) = fRef(end); 
        % Current
        MinRadNorm(dir) = fCurr(1);
        MaxRadNorm(dir) = fCurr(end);  
    end 
end
% Join 
MinMaxJoin = cat(2, MinRad, MaxRad);
MinMaxNormJoin = cat(2, MinRadNorm, MaxRadNorm); 

% Loop through Directories
for i = 1:length(directory)
    % Import data including time, number of steps, and radius
    dataplot = importfile0(directory(i) + "data_plot.txt", 1, inf);
    time = dataplot.VarName1;       % Time 
    steps = dataplot.VarName2;      % Number of steps
    radius = dataplot.VarName3;     % Non-Normalized radius 
    
    % Normalize radius by initial radius 
    radNorm = radius./radius(1);

    % Calculate the point at which the swelling radius is 99% to equilibrium
    % value. Must calculate start point after ramping 
    % Simulation - Nonlinear data 
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
    xlabel("t/\tau")
    ylabel("r")
    hold on 
    h1(i) = plot(time, radius, dashType{i}, 'Color', colors{1,i}, 'LineWidth', 3);
    h2(i) = plot(timeValue(i), radValue(i), '*', 'Color', colors{1,i}, 'LineWidth', 1.5, 'MarkerSize', markSize(i));
    % Limits of x and y axis 
    xlim([1E-1 1E5]) 
    newLim = get(gca,'XLim'); 
    newx = [1E-1, 1, 1E1, 1E2, 1E3, 1E4, 1E5];
    set(gca,'XTick', newx);  
        
    minRadVal = min(MinMaxJoin);
    maxRadVal = max(MinMaxJoin);
    ylim([minRadVal - 0.02, maxRadVal + 0.01])
    
    figure(2) 
    grid on
    xlabel("t/\tau")
    ylabel("r/r_{ini}")
    hold on 
    h5(i) = plot(time, radNorm, dashType{i}, 'Color', colors{1,i}, 'LineWidth', 3);
    h6(i)= plot(timeValue(i), radNormValue(i), '*', 'Color', colors{1,i}, 'LineWidth', 1.5, 'MarkerSize', markSize(i));
    % Limits of x and y axis 
    xlim([1E-1 1E5]) 
    newLim = get(gca,'XLim');
    newx = [1E-1, 1, 1E1, 1E2, 1E3, 1E4, 1E5];
    set(gca,'XTick', newx); 
    set(gca,'fontsize',15, 'XScale', 'log')    

    minRadValNorm = min(MinMaxNormJoin);
    maxRadValNorm = max(MinMaxNormJoin);
    ylim([minRadValNorm - 0.02, maxRadValNorm + 0.01])
    
end 

for a = 1:length(cases)
    load(cases(a))
    % Re-normalizing linear data time 
    NuEq = 0.5 - (n/2)*(1/(l0(a)^2*(l0(a)^3-1)) + n/l0(a)^2 - 2*chi(a)/l0(a)^5 )^(-1);
    TimeNorm = time.*(1/n)*(l0(a)^6/(l0(a)^3 - 1))*(3*(1+2*NuEq)/(2*(3+5*NuEq)));

    % Reference and Current 
    fRef = (f+1)*l0(a);
    fCurr = f+1; 

    % Linear data
    % Non-Normalized
    fDiff = fRef(end) - fRef(1);
    f99 = fRef(end) - fDiff*0.01; 
    % Normalized
    fDiffNorm(a) = fCurr(end) - fCurr(1);
    fNorm99(a) = fCurr(end) - fDiffNorm(a)*0.01; 
    
    minArray = abs(fRef-f99);
    findMin = min(minArray); 
    for aa = 1:length(minArray)
        % Level of precision determined by second value (3) 
        if round(minArray(aa),levPrec) == round(findMin,levPrec) 
            fcloseIndex = aa;
            break
        end
    end

    fRadValue(a) = fRef(fcloseIndex); 
    fRadNormValue(a) = fCurr(fcloseIndex);
    
    % Find time value at the particular index
    fTimeValue(a) = TimeNorm(fcloseIndex);
    
    figure(1) 
    hold on 
    h3(a) = plot(TimeNorm, fRef, dashType{a+ConCol}, 'Color', colors{1,a}, 'LineWidth', 3);
    h4(a) = plot(fTimeValue(a), fRadValue(a), '*', 'Color', colors{1,a}, 'LineWidth', 1.5, 'MarkerSize', markSize(a));
    set(gca,'fontsize',15, 'XScale', 'log')

    figure(2) 
    hold on    
    h7(a) = plot(TimeNorm, fCurr, dashType{a+ConCol}, 'Color', colors{1,a}, 'LineWidth', 3);
    h8(a) = plot(fTimeValue(a), fRadNormValue(a), '*', 'Color', colors{1,a}, 'LineWidth', 1.5, 'MarkerSize', markSize(a));
    set(gca,'fontsize',22, 'XScale', 'log')
end

%%
figure(1) 
for ii = 1:types:length(directory)
    avgRad99 = mean(rad99(1,ii:ii+types-1));
    avgTimeVal = mean(timeValue(1,ii:ii+types-1));
    
    xV = avgTimeVal*ones(1,100); 
    yV = linspace(minRadVal - 0.02,maxRadVal + 0.01,100);
    hold on 
    plot(xV, yV, ':', 'Color', colors{1,end}, 'Linewidth', 2)
    avgFNorm99 = fNorm99(ii);
    fAvgTimeVal = fTimeValue(ii);
    fXV = fAvgTimeVal*ones(1,100);
    plot(fXV, yV, ':', 'Color', colors{1,end}, 'Linewidth', 2)
end

legend([h1],...
      {"$\tilde{\gamma} = 0, \, \lambda_b = \mu_b = 0$",...
      "$\tilde{\gamma} = 1$",...
      "$\lambda_b = \mu_b = 10$",...
      "$\lambda_b = \mu_b = 20$"},...
      'Location', 'Best')

% legend([h1],...
%       {"\gamma = 0.00 " + "\lambda_0 = " + legStretch{1},...
%       "\gamma = 0.01 "  + "\lambda_0 = " + legStretch{2},...
%       "\gamma = 0.00 "  + "\lambda_0 = " + legStretch{3},...
%       "\gamma = 0.01 "  + "\lambda_0 = " + legStretch{4}},...
%       'Location', 'Best')
     
saveas(gcf, 'radius.eps', 'epsc')
saveas(gcf, "radius.fig")

%% 
figure(2) 
for ii = 1:types:length(directory)
    avgRadNorm99 = mean(radNorm99(1,ii:ii+types-1));
    avgTimeVal = mean(timeValue(1,ii:ii+types-1));
    
    xV = avgTimeVal*ones(1,100); 
    yV = linspace(minRadValNorm - 0.1, maxRadValNorm + 0.1,100);
    hold on     
    plot(xV, yV, ':', 'Color', colors{1,end}, 'Linewidth', 2)
    avgFNorm99 = fNorm99(ii);
    fAvgTimeVal = fTimeValue(ii);
    fXV = fAvgTimeVal*ones(1,100);
    plot(fXV, yV, ':', 'Color', colors{1,end}, 'Linewidth', 2)
end 

legend([h5],...
      {"$\tilde{\gamma} = 0, \, \lambda_b = \mu_b = 0$",...
      "$\tilde{\gamma} = 1$",...
      "$\lambda_b = \mu_b = 10$",...
      "$\lambda_b = \mu_b = 20$"},...
      'Location', 'Best')
  
% legend([h5],...
%       {"\gamma = 0.00 " + "\lambda_0 = " + legStretch{1},...
%       "\gamma = 0.01 "  + "\lambda_0 = " + legStretch{2},...
%       "\gamma = 0.00 "  + "\lambda_0 = " + legStretch{3},...
%       "\gamma = 0.01 "  + "\lambda_0 = " + legStretch{4}},...
%       'Location', 'Best')
  
saveas(gcf,'radius_norm.eps','epsc')
saveas(gcf, "radius_norm.fig")

%%

% legend([h1(1), h1(2), h3(1), h3(2)],...
%       {"\lambda_0 = " + legStretch{1} + ", \chi = " + legChi{1},...
%        "\lambda_0 = " + legStretch{3} + ", \chi = " + legChi{3},...
%        "Linear: " + "\lambda_0 = " + legStretch{2} + ", \chi = " + legChi{2},...
%        "Linear: " + "\lambda_0 = " + legStretch{4} + ", \chi = " + legChi{4}},...
%       'Location', 'Best')

% legend([h1(1), h1(2), h1(3), h1(4)],...
%     {"\gamma = 1, \lambda_0 = " + legStretch{1} +  ", \mu_0 = " + legChem{1},...
%      "\gamma = 1, \lambda_0 = " + legStretch{2} +  ", \mu_0 = " + legChem{2},...
%      "\gamma = 0, \lambda_0 = " + legStretch{3} +  ", \mu_0 = " + legChem{3},...
%      "\gamma = 0, \lambda_0 = " + legStretch{4} +  ", \mu_0 = " + legChem{4}},...
%      'Location', 'Best')

% legend([h1(1), h1(2), h1(3), h1(4), h1(5)],...
%     {"\gamma = 0.0, \lambda_0 = " + legStretch{1} +  ", \mu_0 = " + legChem{1},...
%      "\gamma = 0.1, \lambda_0 = " + legStretch{2} +  ", \mu_0 = " + legChem{2},...
%      "\gamma = 0.2, \lambda_0 = " + legStretch{3} +  ", \mu_0 = " + legChem{3},...
%      "\gamma = 0.5, \lambda_0 = " + legStretch{4} +  ", \mu_0 = " + legChem{4},...
%      "\gamma = 1.0, \lambda_0 = " + legStretch{5} +  ", \mu_0 = " + legChem{5}},...
%      'Location', 'Best')
 
% legend([h1(1), h1(2), h1(3), h1(4), h1(5)],...
%         {"\chi = " + legChi{1} +  ", \mu_0 = " + legChem{1},...
%          "\chi = " + legChi{2} +  ", \mu_0 = " + legChem{2},...
%          "\chi = " + legChi{3} +  ", \mu_0 = " + legChem{3},...
%          "\chi = " + legChi{4} +  ", \mu_0 = " + legChem{4},...
%          "\chi = " + legChi{5} +  ", \mu_0 = " + legChem{5}},...
%          'Location', 'Best')
  
% legend([h1(1), h1(2), h1(3), h1(4)],...
%     {"\lambda_0 = " + legStretch{1} +  ", \mu_0 = " + legChem{1},...
%      "\lambda_0 = " + legStretch{2} +  ", \mu_0 = " + legChem{2},...
%      "\lambda_0 = " + legStretch{3} +  ", \mu_0 = " + legChem{3},...
%      "\lambda_0 = " + legStretch{4} +  ", \mu_0 = " + legChem{4}},...
%      'Location', 'Best')

% legend([h1(1), h1(2), h1(3), h1(4)],...
%         {"\gamma = 0.0, \mu_0 = " + legChem{1},...
%          "\gamma = 1.0, \mu_0 = " + legChem{2},...
%          "\lambda_s = \mu_s = 10, \mu_0 = " + legChem{3},...
%          "\lambda_s = \mu_s = 20, \mu_0 = " + legChem{4}},...
%          'Location', 'Best') 

% legend([h1(1), h1(2)],...
%     {"\lambda_0 = " + legStretch{1} +  ", \mu_0 = " + legChem{1},...
%      "\lambda_0 = " + legStretch{2} +  ", \mu_0 = " + legChem{2}},...
%      'Location', 'Best')

% legend([h3(1), h3(2), h3(3), h3(4), h3(5), h3(6)],...
%         {"\lambda_0 = " + legStretch{1} +  ", \mu_0 = " + legChem{1},...
%          "\lambda_0 = " + legStretch{2} +  ", \mu_0 = " + legChem{2},...
%          "\lambda_0 = " + legStretch{3} +  ", \mu_0 = " + legChem{3},...
%          "\lambda_0 = " + legStretch{4} +  ", \mu_0 = " + legChem{4},...
%          "\lambda_0 = " + legStretch{5} +  ", \mu_0 = " + legChem{5},...
%          "\lambda_0 = " + legStretch{6} +  ", \mu_0 = " + legChem{6}},...
%          'Location', 'Best')

% legend([h3(1), h3(2), h3(3), h3(4), h3(5)],...
%         {"\gamma = 0.0, \lambda_0 = " + legStretch{1} +  ", \mu_0 = " + legChem{1},...
%          "\gamma = 0.1, \lambda_0 = " + legStretch{2} +  ", \mu_0 = " + legChem{2},...
%          "\gamma = 0.2, \lambda_0 = " + legStretch{3} +  ", \mu_0 = " + legChem{3},...
%          "\gamma = 0.5, \lambda_0 = " + legStretch{4} +  ", \mu_0 = " + legChem{4},...
%          "\gamma = 1.0, \lambda_0 = " + legStretch{5} +  ", \mu_0 = " + legChem{5}},...
%          'Location', 'Best')

% legend([h3(1), h3(2), h3(3), h3(4), h3(5)],...
%         {"\chi = " + legChi{1} +  ", \mu_0 = " + legChem{1},...
%          "\chi = " + legChi{2} +  ", \mu_0 = " + legChem{2},...
%          "\chi = " + legChi{3} +  ", \mu_0 = " + legChem{3},...
%          "\chi = " + legChi{4} +  ", \mu_0 = " + legChem{4},...
%          "\chi = " + legChi{5} +  ", \mu_0 = " + legChem{5}},...
%          'Location', 'Best')

% legend([h3(1), h3(2), h3(3), h3(4)],...
%     {"\lambda_0 = " + legStretch{1} +  ", \mu_0 = " + legChem{1},...
%      "\lambda_0 = " + legStretch{2} +  ", \mu_0 = " + legChem{2},...
%      "\lambda_0 = " + legStretch{3} +  ", \mu_0 = " + legChem{3},...
%      "\lambda_0 = " + legStretch{4} +  ", \mu_0 = " + legChem{4}},...
%      'Location', 'Best')
 
% legend([h3(1), h3(2), h3(3), h3(4)],...
%         {"\gamma = 0.0, \mu_0 = " + legChem{1},...
%          "\gamma = 1.0, \mu_0 = " + legChem{2},...
%          "\lambda_s = \mu_s = 10, \mu_0 = " + legChem{3},...
%          "\lambda_s = \mu_s = 20, \mu_0 = " + legChem{4}},...
%          'Location', 'Best')  
% 
% legend([h3(1), h3(2)],...
%     {"\gamma = 1, \lambda_0 = " + legStretch{1} +  ", \mu_0 = " + legChem{1},...
%      "\gamma = 1, \lambda_0 = " + legStretch{2} +  ", \mu_0 = " + legChem{2}},...
%      'Location', 'Best')