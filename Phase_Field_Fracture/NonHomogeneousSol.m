clear all
syms c(x) sig l0 c_hom E Gc 
ode = ((2*l0*sig^2)/(c^4*E*Gc) + 1)*c - 4*l0^2*diff(c,x,2) == 1;
cond1 = c(0) == 0;
cond2 = c(Inf) == 1/(2*l0*sig^2/(Gc*c_hom^4*E)+ 1);

conds = [cond1 cond2];
cSol(x) = dsolve(ode,conds)
cSol = simplify(cSol)


%%
syms a b c(x)
f = diff(-a/(c^2) + c^2/2 - 2*b*diff(c,x)^2 -c,x);
Fc = int(f,[x,inf])
Fc = int(f,[-inf,x])

diff(c(x), x, x)