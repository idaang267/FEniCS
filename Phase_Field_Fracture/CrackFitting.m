%%
Delta = [0, .548, .549, .55]'; 
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
a = 0.3006;
b = 0.00189;
c = -0.3006;
d = -0.1254;

x = linspace(0, 350, 350);
y = a*exp(b*x) + c*exp(d*x);

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
