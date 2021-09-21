clear all
close all

% Initial parameters in order to choose our value of delta to aim for in
h = 0.5;    % total strip height is 2*h
mu = 1;     % Shear modulus
Gc = 1;     % Critical fracture energy 

% Calculate the effective Gc 
ell_multi = 5; 
hsize = 0.005; 
ell = ell_multi*hsize;
Gc_e = Gc*(1+ (3/2)*hsize/ell);

syms lambda_a
eq1 = Gc == h*mu*(lambda_a - 1/lambda_a)^2;
eq2 = Gc_e == h*mu*(lambda_a - 1/lambda_a)^2;

sim_lambda = double(solve(eq1, lambda_a));
sim_lambda_e = double(solve(eq2, lambda_a));

sim_delta = h*(sim_lambda(4) - 1)
sim_delta_e = h*(sim_lambda_e(4) - 1)