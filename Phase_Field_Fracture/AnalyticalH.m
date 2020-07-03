clear all 
close all

% Figure 3B replication 

% Theta (rad) 
th = linspace(-pi, pi, 101);

% Fig 3B
n = [1.7, 2, 3];

for ii = 1:length(n)
    % Kappa 
    k = 1 - 1/n(ii); 
    
    % omega [dimensionless function]
    w = sqrt(1 - (k*sin(th)).^2);

    % m
    m = 1 - 1/(2*n(ii));
    
    % Define cos(xi)
    cosxi = 1/(n(ii)*sqrt(2))*(sqrt(1+k*sin(th).^2 - w.*cos(th))./(w+k.*cos(th)));
    % Find xi 
    xi = acos(cosxi);        
    x = cos(xi).^2;
    % Define hypergeometric function F for anglular function H
    for i = 1:length(x)
        F(i) = hypergeom([1/2 - 1/m, 1/2], 3/2 - 1/m, x(i));
    end
    % Angular Function H(th, n)
    H = (-n(ii)^(5/2)/(m^2)).*(w + k*cos(th)).^(2-m).*((m/(2-m))*F - k*sin(xi)); 
    
    % Endpoints
    gamma1 = gamma((2*n(ii)-3)/(2*(2*n(ii)-1)));
    gamma2 = gamma(-1/(2*n(ii)-1));
    H_mid = -((4*n(ii)^(9/2))/(2*n(ii)-1)^2)*(2-1/n(ii))^(1+1/(2*n(ii)))*(1/(n(ii)*(2*n(ii)+1)));
    H_st_en = -((4*n(ii)^(9/2))/(2*n(ii)-1)^2)*n(ii)^(-1-1/(2*n(ii)))*((2*n(ii)-1)/(2*n(ii)+1))*sqrt(pi)*gamma1/gamma2;
    
    X = [th(1), th(51), th(end)];    
    Y = [H_st_en, H_mid, H_st_en];
    clear H_mid H_st_en
    
    figure(1) 
    h1(ii) = plot(th, H, 'LineWidth', 2.5);
    hold on
    scatter(X, Y, 'filled')
    hold on 
    clear X Y
    
end 

%% Plotting
figure(1) 
grid on 
xlabel("$\theta$", 'Interpreter', 'LaTeX')
ylabel("$H(\theta, n)$", 'Interpreter', 'LaTeX')
ax = gca; 
ax.XTick = -pi:pi/2:pi;
ax.XTickLabel = {'-\pi','-\pi/2','0','\pi/2','\pi'};
ax.FontSize = 15; 
xlim([-pi pi])
leg = legend([h1], ["n = 1.7", "n = 2", "n = 3"], 'Location', 'Best');
saveas(gcf, 'Fig3B.eps', 'epsc')
