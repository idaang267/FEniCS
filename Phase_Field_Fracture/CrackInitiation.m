clear all
close all

% Initial parameters in order to choose our value of delta to aim for in
h = 0.25;    % total strip height is 2*h
meshsizeY = 25/2;
mu = 1;     % Shear modulus
Gc = 1;     % Critical fracture energy 
ell_multi = 5; 

% Calculate the effective Gc 
hsizeY = h/meshsizeY; 
ell = ell_multi*hsizeY;
Gc_e = Gc*(1+ (3/8)*hsizeY/ell);

syms lambda_a
eq2 = Gc == h*mu*(lambda_a - 1/lambda_a)^2;
sim_lambda = double(solve(eq2, lambda_a));

sim_delta = h*(sim_lambda(4) - 1)

%%
h = 0.5/2;
mu = 1;
lambda_a = 1 + 0.48/h;
J = h*mu*(lambda_a - 1/lambda_a)^2

%%
Delta = [0, 0.01, .548, .549, .55]'; 
n_steps = [0, 1, 197, 198, 200]'; 

coefs = fit(n_steps, Delta, 'exp2')

a = coefs.a;
b = coefs.b;
c = coefs.c;
d = coefs.d;

x = linspace(0, 200, 200);
y = a*exp(b*x) + c*exp(d*x);

hold on
plot(n_steps, Delta, '*')
plot(x, y, '*')

%%
Dis = TaylorHoodenergies.e03;
ElaEnergy = TaylorHoodenergies.e1;
DisEnergy = TaylorHoodenergies.e00;
TotEnergy = TaylorHoodenergies.e2;

hold on 
plot(Dis(1:36), ElaEnergy(1:36), '*')
plot(Dis(1:36), DisEnergy(1:36), 'LineWidth', 2)
plot(Dis(1:36), TotEnergy(1:36), 'LineWidth', 2)
legend("Elastic","Dissipated","Total")
xlabel("Displacement")
ylabel("Energy")