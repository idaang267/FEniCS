r = 1;
th = linspace(0,2*pi);
Gc = 2.4e6; 
nu = 0.49995;
mu = 1; 
kappa = 2*(1+nu)*mu/(3*(1-2*nu));
E = 9*kappa*mu/(3*kappa+mu);
K_IC = sqrt(Gc*E/(1-nu^2));

u1 = K_IC/mu*sqrt(r/(2*pi))*cos(th/2).*(4/3-cos(th/2).^2);
u2 = K_IC/mu*sqrt(r/(2*pi))*sin(th/2).*(4/3-cos(th/2).^2);

figure
polarplot(th,u1)
hold on 
polarplot(th,u2)


%%
x = linspace(0,1,10);
y = linspace(0,1,10);

Gc = 2.4e6; 
nu = 0.49995;
mu = 1; 
kappa = 2*(1+nu)*mu/(3*(1-2*nu));
E = 9*kappa*mu/(3*kappa+mu);

K_IC = sqrt(Gc*E/(1-nu^2));

r = sqrt(X.^2 + Y.^2);
coshalfx = sqrt((1+X/r)/2);
sinhalfx = sqrt((1-X/r)/2);
com = K_IC/mu*sqrt(r/(2*pi))*(4/3-coshalfx);
u1x = com*coshalfx;
u2x = com*sinhalfx;

surf(X,Y,real(u1x))
surf(X,Y,real(u2x))