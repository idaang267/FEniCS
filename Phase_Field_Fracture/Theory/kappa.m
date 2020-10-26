
% If nu = 0.49995
syms mu
nu = 0.49995;
kppa = 1;
eq = kppa == (2*(1+nu)*mu)/(3*(1-2*nu));
vpasolve(eq, mu)

%%
mu = ans;

mu = 1 
nu = 0.49995;
lmbda = 2*mu*nu/(1-2*nu);
kppa = (2*(1+nu)*mu)/(3*(1-2*nu))
