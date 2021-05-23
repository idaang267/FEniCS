clear all
close all
% Sphere Radius and Plug Radius 
R = 7.1; Rp = 2.5; 

% Angle in relation to the plug 
alpha = Rp/R;
beta = 0;
gamma = pi/2;

points = 4;
for i = 1:points    
    t = (2*pi/points)*i;
    x(i) = sin(alpha)*cos(beta)*cos(gamma)*cos(t)+sin(alpha)*sin(gamma)*sin(t) - cos(alpha)*sin(beta)*cos(gamma);
    y(i) = -sin(alpha)*cos(beta)*sin(gamma)*cos(t)+sin(alpha)*cos(gamma)*sin(t) + cos(alpha)*sin(beta)*sin(gamma);
    z(i) = sin(alpha)*sin(beta)*cos(t)+ cos(alpha)*cos(beta);
end

x_s = [0,R,0,-R,0,0,0];
y_s = [0,0,0,0,0,R,-R];
z_s = [0,0,-R,0,R,0,0];

plot3(R*x,R*y,R*z, 'o')
hold on 
plot3(x_s,y_s,z_s,'o')

%%
[X,Y,Z] = sphere;
hold on 
X2 = X*R;
Y2 = Y*R;
Z2 = Z*R;
surf(X2,Y2,Z2)

  
% 
h = R*(1-cos(th/2));
c = (R^2 - (R-h)^2)^(1/2);


