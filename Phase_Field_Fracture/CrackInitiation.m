% Initial parameters in order to choose our value of delta to aim for in
h = 0.5;    % total strip height is 2*h
mu = 1;     % Shear modulus
Gc = 1;     % Critical fracture energy 

% Unknown stretch 
syms lambda_a
eq = Gc == h*mu*(lambda_a - 1/lambda_a)^2;
th_lambda = double(solve(eq, lambda_a));

% Unknown delta to prescribe 
th_delta = h*(th_lambda(4) - 1)

% Calculate the effective Gc 
meshsizeY = 25;
hsizeY = h/meshsizeY; 
ell_multi = 5; 
ell = ell_multi*hsizeY;
Gc_e = 1;
Gc_sim = Gc_e/((1+ (3/8)*hsizeY/ell));

syms lambda_a
eq2 = Gc_sim == h*mu*(lambda_a - 1/lambda_a)^2;
sim_lambda = double(solve(eq2, lambda_a));

sim_delta = h*(sim_lambda(4) - 1)

%%
Delta = [0, .518, .519, .52]'; 
n_steps = [0, 48, 49, 50]'; 

coefs = fit(n_steps, Delta, 'exp2')

a = 0.4772;
b = 0.001927;
c = -0.4722;
d = -0.7501;

x = linspace(0, 50, 50);
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