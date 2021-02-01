
l0 = 0.01;
E = 10;
G_c = 10;
k = l0*E/(2*G_c);

c_hom = linspace(0,1,100);

eps = sqrt((1 - c_hom)./(k.*c_hom));

stress = E*c_hom.^2.*eps;

plot(eps, c_hom)

plot(eps,stress)

%%
syms c_hom(eps_hom)
l0 = 0.1;
c_hom(eps_hom) = 1/(1 + 2*l0*eps_hom^2);
fplot(c_hom(eps_hom))
