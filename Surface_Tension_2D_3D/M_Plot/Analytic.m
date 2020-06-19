%% 
clear all 

n = 10^-3;
chi = 0.2;
gamma = 1.0;
l0 = 2.0;
J0 = l0^3;
mu0 = n*(1/l0-1/J0) + 1/J0 + log(1 - 1/J0) + chi/(J0^2) %+ n*gamma*2/l0

% No surface tension at the same chemical potential 
syms ls 
Js = ls^3;
Eq1 = (mu0 == n*(1/ls-1/Js) + 1/Js + log(1 - 1/Js) + chi/(Js^2));
l_ans = double(vpasolve(Eq1, ls, 2.0))

%%
n = 10^-3;
chi = 0.4;
gamma = 0.0;
l0 = 2.0;
J0 = l0^3;
mu0 = n*(1/l0-1/J0) + 1/J0 + log(1 - 1/J0) + chi/(J0^2) + n*gamma*2/l0

% No surface tension at the same chemical potential 
syms ls 
Js = ls^3;
Eq1 = (mu0 == n*(1/ls-1/Js) + 1/Js + log(1 - 1/Js) + chi/(Js^2));
l_ans = double(vpasolve(Eq1, ls, 2.0))

%%
clear all
close all
format long

Time = [1E-6, 1E-3, 1E-2, 0.1, 1.0];
n = 10^-3;      % N Omega
chi = [0.4];
l0 = [2.6];
J0 = l0.^3;
NuEq = 0.5 - n/2*(1/(l0^2*(l0^3-1)) + n/l0^2 - 2*chi/l0^5 )^(-1);
TimeNorm = Time.*(1/n)*(l0^6/(l0^3 - 1))*(3*(1+2*NuEq)/(2*(3+5*NuEq)));

% Values we want to fit 
val = [0, 0, TimeNorm(1), TimeNorm(2), TimeNorm(3), TimeNorm(4), TimeNorm(5)]';
% x has to be same length as val
steps = 80;
x = linspace(1,steps,7)';
[curve,gof] = fit(x,val,'exp1')
% from curve we know a and b 
a = curve.a;
b = curve.b;
x_points = linspace(0,steps,steps)';
y = a*exp(b*x_points);
figure
hold on 
scatter(x, val)
plot(x_points, y, '-.')
hold off
%%
steps = 75;
start = 10^-5; 
time(1) = 0;
for i = 2:steps
    time(i) = start*1.6;
    start = time(i);
end 

x = linspace(1,steps,steps)';
[curve,gof] = fit(x,time','exp1')
figure(2)
hold on
plot(x, time, "*")
plot(curve) 


