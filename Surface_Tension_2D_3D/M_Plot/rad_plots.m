clear all
% Comparison of chi = 0.4 with and without surface tension
%directory = ["ex_chi_0.4_l0_1.50/", "ex_chi_0.4_l0_2.00/","ex_chi_0.4_l0_2.50/"];
%directory = ["ex_wo_chi_0.4_l0_1.50/", "ex_wo_chi_0.4_l0_2.00/", "ex_wo_chi_0.4_l0_2.50/"];

% Comparison of varying chi with and without surface tension
%directory = ["ex_chi_0.2_l0_2.00/", "ex_chi_0.4_l0_2.00/"];
%directory = ["ex_wo_chi_0.2_l0_2.00/", "ex_wo_chi_0.4_l0_2.00/"];

%directory = ["exp_chi_0.2_l0_2.00_exp_1.6/", "exp_chi_0.2_l0_2.50_exp_1.6/", "exp_chi_0.2_l0_3.00_exp_1.6/"];
%directory = ["exp_wo_chi_0.2_l0_2.00_exp_1.6/", "exp_wo_chi_0.2_l0_2.50_exp_1.6/", "exp_wo_chi_0.2_l0_3.00_exp_1.6/"];

%directory = ["exp_chi_0.4_l0_2.00_exp_1.6/", "exp_chi_0.4_l0_2.50_exp_1.6/", "exp_chi_0.4_l0_3.00_exp_1.6/"];
%directory = ["exp_wo_chi_0.4_l0_2.00_exp_1.6/", "exp_wo_chi_0.4_l0_2.50_exp_1.6/", "exp_wo_chi_0.4_l0_3.00_exp_1.6/"];

%directory = ["exp_chi_0.2_l0_2.00_exp_1.6/", "exp_chi_0.4_l0_2.00_exp_1.6/"];
%directory = ["exp_wo_chi_0.2_l0_2.00_exp_1.6/", "exp_wo_chi_0.4_l0_2.00_exp_1.6/"]; 

directory = ["alt_chi_0.2_l0_2.00/"];

% Test of alternate formulation
for i = 1:length(directory)
    dataplot = importfile(directory(i) + "data_plot.txt",1,inf);
    time = dataplot.e00;
    steps = dataplot.e1;
    radius = dataplot.e2;
    output = radius./radius(1);
    
    figure(1)
    set(gca,'fontsize',15, 'XScale', 'log')
    xlabel("Time")
    ylabel("R/R_{ini}")
    plot(time,output,'-','LineWidth',2)
    hold on
%   legend("\lambda_0 = 2.0, \mu_0 = -403 \times 10^{-5}", "\lambda_0 = 2.5, \mu_0 = -185 \times 10^{-6}", "\lambda_0 = 3.0, \mu_0 = 534 \times 10^{-6}")
%   legend("\lambda_0 = 2.0, \mu_0 = -906 \times 10^{-6}", "\lambda_0 = 2.5, \mu_0 = 635 \times 10^{-6}", "\lambda_0 = 3.0, \mu_0 = 808 \times 10^{-6}")
%     legend("\gamma = 1.0, \lambda_0 = 2.0, \mu_0 = -403 \times 10^{-5}", ...
%           "\gamma = 1.0, \lambda_0 = 2.5, \mu_0 = -185 \times 10^{-6}", ...
%           "\gamma = 1.0, \lambda_0 = 3.0, \mu_0 = 534 \times 10^{-6}", ...
%           "\gamma = 0.0, \lambda_0 = 2.0, \mu_0 = -503 \times 10^{-5}", ...
%           "\gamma = 0.0, \lambda_0 = 2.5, \mu_0 = -985 \times 10^{-6}", ...
%           "\gamma = 0.0, \lambda_0 = 3.0, \mu_0 = -133 \times 10^{-6}")
%    legend("\gamma = 1.0, \lambda_0 = 2.0, \mu_0 = -403 \times 10^{-5}", ...
%          "\gamma = 1.0, \lambda_0 = 2.5, \mu_0 = -185 \times 10^{-6}", ...
%          "\gamma = 1.0, \lambda_0 = 3.0, \mu_0 = 534 \times 10^{-6}", ...
%          "\gamma = 0.0, \lambda_0 = 2.0, \mu_0 = -191 \times 10^{-5}", ...
%          "\gamma = 0.0, \lambda_0 = 2.5, \mu_0 = -165 \times 10^{-6}", ...
%          "\gamma = 0.0, \lambda_0 = 3.0, \mu_0 = 142 \times 10^{-6}")
% 
%     legend("\gamma = 1.0, \chi = 0.2, \mu_0 = -403 \times 10^{-5}",...
%            "\gamma = 1.0, \chi = 0.4, \mu_0 = -906 \times 10^{-6}",...
%            "\gamma = 0.0, \chi = 0.2, \mu_0 = -503 \times 10^{-5}",...
%            "\gamma = 0.0, \chi = 0.4, \mu_0 = -191 \times 10^{-5}")
%     
    saveas(gcf,'radius_norm.png')

    figure(2) 
    set(gca,'fontsize',15, 'XScale', 'log')
    xlabel("Time")
    ylabel("Radius")
    plot(time,radius,'-','LineWidth',2)
    hold on       
%    legend("\lambda_0 = 2.0, \mu_0 = -403 \times 10^{-5}", "\lambda_0 = 2.5, \mu_0 = -185 \times 10^{-6}", "\lambda_0 = 3.0, \mu_0 = 534 \times 10^{-6}")
%    legend("\lambda_0 = 2.0, \mu_0 = -906 \times 10^{-6}", "\lambda_0 = 2.5, \mu_0 = 635 \times 10^{-6}", "\lambda_0 = 3.0, \mu_0 = 808 \times 10^{-6}")
%     legend("\gamma = 1.0, \lambda_0 = 2.0, \mu_0 = -403 \times 10^{-5}", ...
%           "\gamma = 1.0, \lambda_0 = 2.5, \mu_0 = -185 \times 10^{-6}", ...
%           "\gamma = 1.0, \lambda_0 = 3.0, \mu_0 = 534 \times 10^{-6}", ...
%           "\gamma = 0.0, \lambda_0 = 2.0, \mu_0 = -503 \times 10^{-5}", ...
%           "\gamma = 0.0, \lambda_0 = 2.5, \mu_0 = -985 \times 10^{-6}", ...
%           "\gamma = 0.0, \lambda_0 = 3.0, \mu_0 = -133 \times 10^{-6}")
%     legend("\gamma = 1.0, \lambda_0 = 2.0, \mu_0 = -403 \times 10^{-5}", ...
%          "\gamma = 1.0, \lambda_0 = 2.5, \mu_0 = -185 \times 10^{-6}", ...
%          "\gamma = 1.0, \lambda_0 = 3.0, \mu_0 = 534 \times 10^{-6}", ...
%          "\gamma = 0.0, \lambda_0 = 2.0, \mu_0 = -191 \times 10^{-5}", ...
%          "\gamma = 0.0, \lambda_0 = 2.5, \mu_0 = -165 \times 10^{-6}", ...
%          "\gamma = 0.0, \lambda_0 = 3.0, \mu_0 = 142 \times 10^{-6}")
%       legend("\gamma = 1.0, \chi = 0.2, \mu_0 = -403 \times 10^{-5}",...
%            "\gamma = 1.0, \chi = 0.4, \mu_0 = -906 \times 10^{-6}",...
%            "\gamma = 0.0, \chi = 0.2, \mu_0 = -503 \times 10^{-5}",...
%            "\gamma = 0.0, \chi = 0.4, \mu_0 = -191 \times 10^{-5}")
%       
       saveas(gcf,'radius.png')

end 

%%
n = 10^-3;
gamma = 0.0;
chi = 0.4;
l0 = 3.0;
J0 = l0^3;
mu0 = n*(1/l0-1/J0) + 1/J0 + log(1 - 1/J0) + chi/(J0^2)+ n*gamma*2/l0

