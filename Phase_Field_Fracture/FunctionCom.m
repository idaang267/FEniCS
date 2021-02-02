close all

% Alpha only goes from 0 to 1 
alpha = linspace(0,1,200);

a = (1-alpha).^2;
b = (1-alpha).^3;

%% Test
q = 1024;
m = 6;
mod = (a./(a + m*alpha.*(1+q*alpha))).^2;
m = 6;
mod2 = (a./(a + m*alpha.*(1+q*alpha))).^2;
m = 24;
mod3 = (a./(a + m*alpha.*(1+q*alpha))).^2;

g = 1.1;
gd = g*(1-((g-1)/g).^(alpha.^2));
gm = 1 - 3*alpha.^2+2*alpha.^3;

hold on 
h0 = plot(alpha, a, 'Linewidth', 2.5);
h1 = plot(alpha, b, 'Linewidth', 2.5);
h2 = plot(alpha, gd, 'Linewidth', 2.5); 
h3 = plot(alpha, gm, 'Linewidth', 2.5); 
h4 = plot(alpha, mod3, 'Linewidth', 2.5); 

% h2 = plot(alpha, gd, 'Linewidth', 2.5); 
% h3 = plot(alpha, gm, 'Linewidth', 2.5); 
leg = legend([h0, h1, h2, h3, h4],["a", "b", "gd", "gm", "q=64, m=12"], 'Location', 'Best');

%% Macauley Bracket and Heaviside Function
mac = (alpha + abs(alpha))/2;
Heav = (alpha + abs(alpha))./(2.*alpha);

% Shifted Heaviside
k_ell = .99;
Heav_shift = (alpha-k_ell + abs(alpha-k_ell))./(2.*(alpha-k_ell));

% Logistic Function
log_fcn = 1./(1+exp(-200*(alpha-0.9)));
% Hyperbolic Tangent 
htan_fcn = (exp(alpha) - exp(-alpha))./(exp(alpha) + exp(-alpha));
% Gen
gen = (alpha)./sqrt(1 + alpha.^2);

% mu = 0.94*ones([1,200]);
% beta = 0.01; 
% t_exp = -(alpha-mu)/beta;
% den = exp(t_exp) + 1;
% P = 1./den;

% h00 = plot(alpha, a, 'Linewidth', 2.5);
hold on
h0 = plot(alpha, b, 'Linewidth', 2.5);
% h1 = plot(alpha, mac, 'Linewidth', 2.5); 
% h2 = plot(alpha, mac_shift,'Linewidth', 2.5);
% h3 = plot(alpha, Heav,'Linewidth', 2.5);
% h1 = plot(alpha, gen, 'Linewidth', 2.5); 
hold on 
h2 = plot(alpha, gd, 'Linewidth', 2.5); 
h3 = plot(alpha, gm, 'Linewidth', 2.5); 

% h3 = plot(alpha, htan_fcn, 'Linewidth', 2.5); 
h4 = plot(alpha, Heav_shift,'Linewidth', 2.5);
% h5 = plot(alpha, P,'Linewidth', 2.5);

leg = legend([h0, h2, h4],["b", "Logistic", "Heaviside Shift"], 'Location', 'Best');