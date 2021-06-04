%%
Delta = [0, .528, .529, .53]'; 
n_steps = [0, 298, 299, 300]'; 

coefs = fit(n_steps, Delta, 'exp2')

a = coefs.a;
b = coefs.b;
c = coefs.c;
d = coefs.d;

x = linspace(0, 300, 300);
y = a*exp(b*x) + c*exp(d*x);

hold on
plot(n_steps, Delta, '*')
plot(x, y, '*')

%%
Delta = [0, 0.1, 2.58, 2.59, 2.6]'; 
n_steps = [0, 1, 198 , 199, 200]'; 

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
