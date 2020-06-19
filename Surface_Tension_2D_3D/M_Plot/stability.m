%% 
close all
clear all
%% 

lambda_0_01 = [3.0, 2.0, 1.9, 1.8, 1.7, 1.65, 1.6, 1.55, 1.5, 1.45, 1.4, 1.35];
gamma_01 = [0, 0, 0, 0, 0, 0.2, 0.2, 1, 5, 6.9, 13.7, 19]; 

for i = 1:length(lambda_0_01)
    plot(gamma_01(i), lambda_0_01(i), '*', 'MarkerSize', 10)
    hold on 
end

%% 

% Fit for exponential with two terms 
fo1 = fitoptions('Method','NonlinearLeastSquares',...
               'Robust','Bisquare',...
               'Algorithm', 'Trust-Region',...
               'MaxIter', 1000);
% For chi = 0.1
lambda_0_01 = [3.0, 2.0, 1.9, 1.8, 1.7, 1.65, 1.6, 1.55, 1.5, 1.45, 1.4, 1.35]';
gamma_01 = [0, 0, 0, 0, 0, 0.2, 0.2, 1, 5, 6.9, 13.7, 19]'; 
[curve1,gof1] = fit(gamma_01, lambda_0_01,'exp2', fo1)
a = 0.2985; 
b = -6.842;
c = 1.551;
d = -0.007551;
x_gamma_01 = linspace(10^-3,100,1E5); 
y_01 = a*exp(b*x_gamma_01) + c*exp(d*x_gamma_01);

% Set R Square value 
legend_01 = num2str(gof1.rsquare, '%.4f\n');

% For chi = 0.2
fo2 = fitoptions('Method','NonlinearLeastSquares',...
               'Robust','Bisquare',...
               'Algorithm', 'Trust-Region',...
               'MaxIter', 500,...
                'StartPoint', [a, b, c, d]);
            
lambda_0_02 = [2.5, 2.0, 1.9, 1.8, 1.7, 1.6, 1.55, 1.5, 1.45, 1.4, 1.35]';
gamma_02 = [0, 0, 0, 0, 0, 0.1, 0.5, 1, 4.9, 6, 12]'; 
[curve2,gof2] = fit(gamma_02, lambda_0_02, 'exp2', fo2)
a = 0.3186; 
b = -15.06;
c = 1.531;
d = -0.01144;
x_gamma_02 = linspace(10^-3,100,1E5); 
y_02 = a*exp(b*x_gamma_02) + c*exp(d*x_gamma_02);
        
legend_02 = num2str(gof2.rsquare, '%.4f\n');

% For chi = 0.3
fo3 = fitoptions('Method','NonlinearLeastSquares',...
                 'Robust','LAR',...
                 'Algorithm', 'Levenberg-Marquardt',...
                 'MaxIter', 100,...
                 'StartPoint', [0.329, -18, 1.3841, -0.0227]);
lambda_0_03 = [3.0, 2.0, 1.9, 1.8, 1.7, 1.6, 1.55, 1.5, 1.45, 1.4, 1.35, 1.3]';
gamma_03 = [ 0, 0, 0, 0, 0, 0, 0, 0.5, 0.8, 3, 8, 18]'; 
[curve3,gof3] = fit(gamma_03, lambda_0_03,'exp2', fo3)
a = 0.3426; 
b = -17.58;
c = 1.457;
d = -0.006349;
x_gamma_03 = linspace(10^-3,100,1E5); 
y_03 = a*exp(b*x_gamma_03) + c*exp(d*x_gamma_03);

legend_03 = num2str(gof3.rsquare, '%.4f\n');

% For chi = 0.4
fo4 = fitoptions('Method','NonlinearLeastSquares',...
               'Robust','LAR',...
               'Algorithm', 'Levenberg-Marquardt',...
               'MaxIter', 300,...
                'StartPoint', [a, b, c, d]);
lambda_0_04 = [3.0, 2.0, 1.9, 1.8, 1.7, 1.6, 1.5, 1.45, 1.4, 1.35, 1.3]';
gamma_04 = [0, 0, 0, 0, 0, 0, 0, 0.1, 0.7, 3, 8]'; 
[curve4,gof4] = fit(gamma_04, lambda_0_04,'exp2', fo4)
a = 0.39; 
b = -22.42;
c = 1.41;
d = -0.01015;
x_gamma_04 = linspace(10^-3,100,1E5); 
y_04 = a*exp(b*x_gamma_04) + c*exp(d*x_gamma_04);

legend_04 = num2str(gof4.rsquare, '%.4f\n');

% For chi = 0.5
fo5 = fitoptions('Method','NonlinearLeastSquares',...
               'Robust','LAR',...
               'Algorithm', 'Levenberg-Marquardt',...
               'MaxIter', 300,...
                'StartPoint', [a, -15, c, d]);
lambda_0_05 = [3.0, 2.0, 1.7, 1.6, 1.5, 1.45, 1.4, 1.35, 1.3, 1.28, 1.26, 1.25, 1.24, 1.22, 1.2]';
gamma_05 = [0, 0, 0, 0, 0, 0, 0, 0, 1, 3.9, 8.3, 8.8, 12.5, 17.9, 19.5]'; 
[curve5,gof5] = fit(gamma_05, lambda_0_05,'exp2', fo5)
a = 0.2935; 
b = -15;
c = 1.306;
d = -0.004281;
x_gamma_05 = linspace(10^-3,100,1E5); 
y_05 = a*exp(b*x_gamma_05) + c*exp(d*x_gamma_05);

legend_05 = num2str(gof5.rsquare, '%.4f\n');
        
% For chi = 0.6
fo6 = fitoptions('Method','NonlinearLeastSquares',...
               'Robust','Bisquare',...
               'Algorithm', 'Trust-Region',...
               'MaxIter', 1000,...
                'StartPoint', [a, b, c, d]);       
        
lambda_0_06 = [ 2.0, 1.5, 1.4, 1.3, 1.28, 1.26, 1.25, 1.24, 1.22, 1.2]';
gamma_06 = [0, 0, 0, 0, 0.1, 0.3, 0.2, 0.9, 5.9, 9.9]';
[curve6,gof6] = fit(gamma_06, lambda_0_06,'exp2', fo6)
a = 0.1496; 
b = -16.35;
c = 1.25;
d = -0.004188;
x_gamma_06 = linspace(10^-3,100,1E5); 
y_06 = a*exp(b*x_gamma_06) + c*exp(d*x_gamma_06);

legend_06 = num2str(gof6.rsquare, '%.4f\n'); 


%% 
figure(1)
hold on 
a = area(x_gamma_01, y_01, 'LineStyle', 'none');
b = plot(x_gamma_01, y_01, 'k-', 'LineWidth', 2);
c = plot(gamma_01, lambda_0_01, 'o', 'MarkerSize', 6);
hold off 

xlabel("$\tilde{\gamma}$", 'Interpreter', 'LaTeX')
ylabel("\lambda_0")
leg = legend([b, c], "Fitted Curve: " + legend_01, "Raw Data")
set(gca,'fontsize',20, 'XScale', 'log')
ylim([1 2])
xlim([10^-3 100]) 
newx = [1E-3, 1E-2, 1E-1, 1, 1E1, 1E2];
set(gca,'XTick', newx); 
saveas(gcf, 'fit_1.eps', 'epsc')
saveas(gcf, "fit_1.fig")

figure(2) 
hold on 
plot(x_gamma_02, y_02, '-', 'LineWidth', 2);
plot(gamma_02, lambda_0_02, '*', 'MarkerSize', 6, 'LineWidth', 2);
legend("Fitted Curve: " + legend_02, "Raw Data")
set(gca,'fontsize',20, 'XScale', 'log')
xlim([10^-3 10^2]) 
hold off
      
figure(3)
hold on 
plot(x_gamma_03, y_03, '-', 'LineWidth', 2);
plot(gamma_03, lambda_0_03, '*', 'MarkerSize', 6, 'LineWidth', 2);
legend("Fitted Curve: " + legend_03, "Raw Data")
set(gca,'fontsize',20, 'XScale', 'log')
hold off

figure(4) 
hold on 
plot(x_gamma_04, y_04, '-', 'LineWidth', 2);
plot(gamma_04, lambda_0_04, '*', 'MarkerSize', 6, 'LineWidth', 2);
legend("Fitted Curve: " + legend_04, "Raw Data")
set(gca,'fontsize',20, 'XScale', 'log')
hold off

figure(5) 
hold on 
a = area(x_gamma_05, y_05, 'LineStyle', 'none');
b = plot(x_gamma_05, y_05, 'k-', 'LineWidth', 2);
c = plot(gamma_05, lambda_0_05, '*', 'MarkerSize', 6, 'LineWidth', 2);
hold off
minY = min(y_05) - 0.1;
maxY = max(y_05) + 0.1; 
ylim([minY maxY])
xlim([10^-3 10^1]) 
xlabel("$\tilde{\gamma}$", 'Interpreter', 'LaTeX')
ylabel("\lambda_0")
leg5 = legend([b, c],"Fitted Curve: " + legend_05, "Raw Data")
set(gca,'fontsize',20, 'XScale', 'log')

% saveas(gcf, 'fit_5.eps', 'epsc')
% saveas(gcf, "fit_5.fig")

figure(6) 
hold on
a = area(x_gamma_06, y_06, 'LineStyle', 'none');
b = plot(x_gamma_06, y_06, 'k-', 'LineWidth', 2);
c = plot(gamma_06, lambda_0_06, '*', 'MarkerSize', 6, 'LineWidth', 2);
hold off
minY = min(y_06) - 0.1;
maxY = max(y_06) + 0.1; 
ylim([minY maxY])
xlim([10^-3 10^1]) 
xlabel("$\tilde{\gamma}$", 'Interpreter', 'LaTeX')
ylabel("\lambda_0")
leg = legend([b, c],"Fitted Curve: " + legend_06, "Raw Data")
set(gca,'fontsize',20, 'XScale', 'log')

% saveas(gcf, 'fit_6.eps', 'epsc')
% saveas(gcf, "fit_6.fig")

%% 
figure(8) 
grid
hold on
plot(x_gamma_01, y_01, 'k-', 'LineWidth', 2);
plot(x_gamma_02, y_02, '-', 'LineWidth', 2);
plot(x_gamma_03, y_03, '-', 'LineWidth', 2);
plot(x_gamma_04, y_04, '-', 'LineWidth', 2);
plot(x_gamma_05, y_05, '-', 'LineWidth', 2);
plot(x_gamma_06, y_06, '-', 'LineWidth', 2);
hold off
xlim([1E-3 25])
xlabel("$\tilde{\gamma}$", 'Interpreter', 'LaTeX')
ylabel("\lambda_0")
leg = legend("0.1", "0.2", "0.3", "0.4","0.5", "0.6",...
             'Location', 'Best');
title(leg,'\chi')
set(gca,'fontsize',20, 'XScale', 'log')
newx = [1E-3, 1E-2, 1E-1, 1, 1E1, 1E2];
set(gca,'XTick', newx); 
saveas(gcf, 'Combined.eps', 'epsc')
saveas(gcf, "Combined.fig")

%%
figure
% Gamma
X = cat(1,x_gamma_01,x_gamma_02,x_gamma_03,x_gamma_04, x_gamma_05, x_gamma_06);
% Chi 
Y = [0.1*ones(1,1E5); 0.2*ones(1,1E5);...
     0.3*ones(1,1E5); 0.4*ones(1,1E5);...
     0.5*ones(1,1E5); 0.6*ones(1,1E5)];
% Lambda 
Z = cat(1, y_01, y_02, y_03, y_04, y_05, y_06);
surf(X,Y,Z,'FaceAlpha',0.5)
xlim([1E-4 100])
ylim([0 0.7])
zlim([1 2])
% chi_01 = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]';
% plot3(gamma_01, lambda_0_01, chi_01, 'o')
xlabel("$\tilde{\gamma}$", 'Interpreter', 'LaTeX')
ylabel("\chi")
zlabel("\lambda_0")
set(gca,'fontsize',20, 'XScale', 'log')
newx = [1E-4, 1E-2, 1, 1E2];
newy = [0, 0.2, 0.4, 0.6];
set(gca,'XTick', newx); 
set(gca,'YTick', newy); 

saveas(gcf, 'surf.eps', 'epsc')
saveas(gcf, "surf.fig")

%%
figure
plot(gamma_01, lambda_0_01, '-o', 'MarkerSize', 6, 'LineWidth', 2);
hold on 
plot(gamma_02, lambda_0_02, '-o', 'MarkerSize', 6, 'LineWidth', 2);
hold on 
plot(gamma_03, lambda_0_03, '-o', 'MarkerSize', 6, 'LineWidth', 2);
hold on 
plot(gamma_04, lambda_0_04, '-o', 'MarkerSize', 6, 'LineWidth', 2);
hold on 
plot(gamma_05, lambda_0_05, '-o', 'MarkerSize', 6, 'LineWidth', 2);
hold on
plot(gamma_06, lambda_0_06, '-o', 'MarkerSize', 6, 'LineWidth', 2);

leg = legend("0.1","0.2", "0.3", "0.4", "0.5", "0.6",'Location', 'Best');
title(leg,'\chi');
xlabel("\gamma")
ylabel("\lambda_0 ")
set(gca,'fontsize',20)


