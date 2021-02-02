clear all
close all
%%
r = 1;
th = linspace(-pi, pi, 101);
y1 = sqrt(r/(2*pi)).*(7/3*cos(th./2) - cos(3*th./2)) + ...
     (1/(4*sqrt(2*pi)))^2*(-1/15*log(r)-52/45*(log(r)+3/4*sin(th).^2) - 103/48*cos(th) + 26/15*cos(2*th)-3/16*cos(3*th));
y2 = sqrt(r/(2*pi)).*(13/3*sin(th./2) - sin(3*th./2)) + ...
    (1/(4*sqrt(2*pi)))^2*(th/15 - 52/45*(th/4 - 3/8*sin(2*th)) - 61/48*sin(th) + 26/15*sin(2*th) - 3/16*sin(3*th)) ;

figure(1)
hold on 
plot(th, y1, 'linewidth', 2.5)
plot(th, y2, 'linewidth', 2.5)
xline(-pi, '--k')
xline(-3*pi/4, '--k')
xline(-2*pi/4, '--k')
xline(-pi/4, '--k')
xline(0, '--k')
xline(pi/4, '--k')
xline(2*pi/4, '--k')
xline(3*pi/4, '--k')
xline(pi, '--k')

% h0 is the height 
h0 = 100000.0;
% Amount of displacement 
Delta = 100; 
% stretch in direction of displacement 
lambda_a = 1 + Delta/h0;

% Known a for specific geometry 
a = 2*sqrt(h0/pi)*(lambda_a - inv(lambda_a));
c = 1.55;
%%
k_I = 0.001;
mu = 1;
scale = 1.0;

test_par = k_I/(4*mu); 

%% Define LHS section - discontinuous section 
x_l_b = linspace(-0.5, -0.5, 101);
y_l_b = linspace(-0.5, -0, 101);
% Convert coordinates
r_l_b = sqrt(x_l_b.^2 + y_l_b.^2);
th_l_b = atan(y_l_b./x_l_b) - pi;

x_l_t = linspace(-0.5, -0.5, 101);
y_l_t = linspace(0, 0.5, 101);
% Convert coordinates
r_l_t = sqrt(x_l_t.^2 + y_l_t.^2);
th_l_t = atan(y_l_t./x_l_t) + pi;

y1_l_b = test_par*sqrt(r_l_b/(2*pi)).*(7/3*cos(th_l_b./2) - cos(3*th_l_b./2)) + ...
         (test_par/(sqrt(2*pi)))^2*(-1/15*log(r_l_b)-52/45*(log(r_l_b)+3/4*sin(th_l_b).^2) - 103/48*cos(th_l_b) + 26/15*cos(2*th_l_b)-3/16*cos(3*th_l_b));
%c*r.*cos(th-pi) - x_l_b;
y2_l_b = test_par*sqrt(r_l_b/(2*pi)).*(13/3*sin(th_l_b./2) - sin(3*th_l_b./2)) + ...
         (test_par/(sqrt(2*pi)))^2*(th_l_b/15 - 52/45*(th_l_b/4 - 3/8*sin(2*th_l_b)) - 61/48*sin(th_l_b) + 26/15*sin(2*th_l_b) - 3/16*sin(3*th_l_b)) ;
%a*sqrt(r).*sin((th-pi)./2) - y_l_b;
y1_l_t = test_par*sqrt(r_l_t/(2*pi)).*(7/3*cos(th_l_t./2) - cos(3*th_l_t./2)) + ...
         (test_par/(sqrt(2*pi)))^2*(-1/15*log(r_l_t) - 52/45*(log(r_l_t) + 3/4*sin(th_l_t).^2) - 103/48*cos(th_l_t) + 26/15*cos(2*th_l_t) -3/16*cos(3*th_l_t));
% c*r.*cos(th+pi) - x_l_t;
y2_l_t = test_par*sqrt(r_l_t/(2*pi)).*(13/3*sin(th_l_t./2) - sin(3*th_l_t./2)) + ...
         (test_par/(sqrt(2*pi)))^2*(th_l_t/15 - 52/45*(th_l_t/4 - 3/8*sin(2*th_l_t)) - 61/48*sin(th_l_t) + 26/15*sin(2*th_l_t) - 3/16*sin(3*th_l_t)) ;
%a*sqrt(r).*sin((th+pi)./2) - y_l_t;

%% Define Top section - continuous
x_t_l = linspace(-0.5, 0, 101);
y_t_l = linspace(0.5, 0.5, 101);
% Convert coordinates
r_t_l = sqrt(x_t_l.^2 + y_t_l.^2);
th_t_l = atan(y_t_l./x_t_l);
th_t_l(101) = -th_t_l(101);
th_t_l = th_t_l + pi; 

x_t_r = linspace(0, 0.5, 101);
y_t_r = linspace(0.5, 0.5, 101);
% Convert coordinates
r_t_r = sqrt(x_t_r.^2 + y_t_r.^2);
th_t_r = atan(y_t_r./x_t_r);

y1_t_l = test_par*sqrt(r_t_l/(2*pi)).*(7/3*cos(th_t_l./2) - cos(3*th_t_l./2)) + ...
         test_par^2*(1/(sqrt(2*pi)))^2*(-1/15*log(r_t_l) - 52/45*(log(r_t_l) + 3/4*sin(th_t_l).^2) - 103/48*cos(th_t_l) + 26/15*cos(2*th_t_l) -3/16*cos(3*th_t_l));
%flip(c*r.*cos(th + pi) - x_t_l);
y2_t_l = test_par*sqrt(r_t_l/(2*pi)).*(13/3*sin(th_t_l./2) - sin(3*th_t_l./2)) + ...
         test_par^2*(1/(sqrt(2*pi)))^2*(th_t_l/15 - 52/45*(th_t_l/4 - 3/8*sin(2*th_t_l)) - 61/48*sin(th_t_l) + 26/15*sin(2*th_t_l) - 3/16*sin(3*th_t_l)) ;
%flip(a*sqrt(r).*sin((th + pi)./2) - y_t_l);
y1_t_r = test_par*sqrt(r_t_r/(2*pi)).*(7/3*cos(th_t_r./2) - cos(3*th_t_r./2)) + ...
         test_par^2*(1/(sqrt(2*pi)))^2*(-1/15*log(r_t_r) - 52/45*(log(r_t_r) + 3/4*sin(th_t_r).^2) - 103/48*cos(th_t_r) + 26/15*cos(2*th_t_r) -3/16*cos(3*th_t_r));
%c*r.*cos(th) - x_t_r;
y2_t_r = test_par*sqrt(r_t_r/(2*pi)).*(13/3*sin(th_t_r./2) - sin(3*th_t_r./2)) + ...
         test_par^2*(1/(sqrt(2*pi)))^2*(th_t_r/15 - 52/45*(th_t_r/4 - 3/8*sin(2*th_t_r)) - 61/48*sin(th_t_r) + 26/15*sin(2*th_t_r) - 3/16*sin(3*th_t_r)) ;
%a*sqrt(r).*sin((th)./2) - y_t_r;

%% Define RHS - continuous
x_rhs = linspace(0.5, 0.5, 101);
y_rhs = linspace(0.5, -0.5, 101);
% Convert coordinates
r_r = sqrt(x_rhs.^2 + y_rhs.^2);
th_r = atan(y_rhs./x_rhs);
thdeg = rad2deg(th_r)

y1_rhs = test_par*sqrt(r_r/(2*pi)).*(7/3*cos(th_r./2) - cos(3*th_r./2)) + ...
         test_par^2*(1/(sqrt(2*pi)))^2*(-1/15*log(r_r) - 52/45*(log(r_r) + 3/4*sin(th_r).^2) - 103/48*cos(th_r) + 26/15*cos(2*th_r) -3/16*cos(3*th_r));
%c*r.*cos(th) - x_rhs;
y2_rhs = test_par*sqrt(r_r/(2*pi)).*(13/3*sin(th_r./2) - sin(3*th_r./2)) + ...
         test_par^2*(1/(sqrt(2*pi)))^2*(th_r/15 - 52/45*(th_r/4 - 3/8*sin(2*th_r)) - 61/48*sin(th_r) + 26/15*sin(2*th_r) - 3/16*sin(3*th_r)) ;
%flip(a*sqrt(r).*sin((th + pi)./2) - y_t_l);
%a*sqrt(r).*sin((th)./2) - y_rhs;

%% Define Bottom - continuous 
x_b_r = linspace(0.5, 0, 101);
y_b_r = linspace(-0.5, -0.5, 101);
% Convert coordinates
r_b_r = sqrt(x_b_r.^2 + y_b_r.^2);
th_b_r = atan(y_b_r./x_b_r);

x_b_l = linspace(0, -0.5, 101);
y_b_l = linspace(-0.5, -0.5, 101);
% Convert coordinates
r_b_l = sqrt(x_b_l.^2 + y_b_l.^2);
th_b_l = atan(y_b_l./x_b_l);
th_b_l(1) = -th_b_l(1);
th_b_l = th_b_l - pi;

y1_b_r = test_par*sqrt(r_b_r/(2*pi)).*(7/3*cos(th_b_r./2) - cos(3*th_b_r./2)) + ...
         test_par^2*(1/(sqrt(2*pi)))^2*(-1/15*log(r_b_r) - 52/45*(log(r_b_r) + 3/4*sin(th_b_r).^2) - 103/48*cos(th_b_r) + 26/15*cos(2*th_b_r) -3/16*cos(3*th_b_r));
y2_b_r = test_par*sqrt(r_b_r/(2*pi)).*(13/3*sin(th_b_r./2) - sin(3*th_b_r./2)) + ...
         test_par^2*(1/(sqrt(2*pi)))^2*(th_b_r/15 - 52/45*(th_b_r/4 - 3/8*sin(2*th_b_r)) - 61/48*sin(th_b_r) + 26/15*sin(2*th_b_r) - 3/16*sin(3*th_b_r)) ;

y1_b_l = test_par*sqrt(r_b_l/(2*pi)).*(7/3*cos(th_b_l./2) - cos(3*th_b_l./2)) + ...
         test_par^2*(1/(sqrt(2*pi)))^2*(-1/15*log(r_b_l) - 52/45*(log(r_b_l) + 3/4*sin(th_b_l).^2) - 103/48*cos(th_b_l) + 26/15*cos(2*th_b_l) -3/16*cos(3*th_b_l));
y2_b_l = test_par*sqrt(r_b_l/(2*pi)).*(13/3*sin(th_b_l./2) - sin(3*th_b_l./2)) + ...
         test_par^2*(1/(sqrt(2*pi)))^2*(th_b_l/15 - 52/45*(th_b_l/4 - 3/8*sin(2*th_b_l)) - 61/48*sin(th_b_l) + 26/15*sin(2*th_b_l) - 3/16*sin(3*th_b_l)) ;

%%
figure(2)
hold on
plot(linspace(-1.5, -1, 101), scale*y1_b_r, 'r', 'LineWidth', 2.5)
plot(linspace(-1.5, -1, 101), scale*y2_b_r, 'b', 'LineWidth', 2.5)
plot(linspace(-1, -0.5, 101), scale*y1_b_l, 'r', 'LineWidth', 2.5)
plot(linspace(-1, -0.5, 101), scale*y2_b_l, 'b', 'LineWidth', 2.5)
xline(-0.5, '--k')
plot(y_l_t, scale*y1_l_t, 'r', 'LineWidth', 2.5)
plot(y_l_t, scale*y2_l_t, 'b', 'LineWidth', 2.5)
h1 = plot(y_l_b, scale*y1_l_b, 'r', 'LineWidth', 2.5);
h2 = plot(y_l_b, scale*y2_l_b, 'b', 'LineWidth', 2.5);
xline(0.5, '--k') % top 
plot(linspace(0.5, 1., 101), scale*y1_t_l, 'r', 'LineWidth', 2.5)
plot(linspace(0.5, 1., 101), scale*y2_t_l, 'b', 'LineWidth', 2.5)
plot(linspace(1.0, 1.5, 101), scale*y1_t_r, 'r', 'LineWidth', 2.5)
plot(linspace(1.0, 1.5, 101), scale*y2_t_r, 'b', 'LineWidth', 2.5)
xline(1.5, '--k') % rhs
plot(linspace(1.5, 2.5, 101), scale*y1_rhs, 'r', 'LineWidth', 2.5)
plot(linspace(1.5, 2.5, 101), scale*y2_rhs, 'b', 'LineWidth', 2.5)

leg = legend([h1, h2],["u1", "u2"], 'Location', 'Best');

%%
x_fea = KI0S1.Points0;
ux_fea = KI0S1.Disp0;
uy_fea = KI0S1.Disp1;
k_I = 0.001;
mu = 1;

par = k_I/(4*mu); 

x = linspace(-0.5,0,201);
y = linspace(0.0,0.0,201);

r = sqrt(x.^2 + y.^2);
th = atan(y./x)+pi;
r = r(1,1:200);
th = th(1,1:200);
ux = par*sqrt(r/(2*pi)).*(7/3*cos(th./2) - cos(3*th./2)) + ...
     par^2*(1/(sqrt(2*pi)))^2*(-1/15*log(r) - 52/45*(log(r) + 3/4*sin(th).^2) - 103/48*cos(th) + 26/15*cos(2*th) -3/16*cos(3*th));
uy = par*sqrt(r/(2*pi)).*(13/3*sin(th./2) - sin(3*th./2)) + ...
     par^2*(1/(sqrt(2*pi)))^2*(th/15 - 52/45*(th/4 - 3/8*sin(2*th)) - 61/48*sin(th) + 26/15*sin(2*th) - 3/16*sin(3*th)) ;
% ux - ux(0, 0) and uy - uy(0, 0)
uy_0 = uy - uy(end);
ux_0 = ux - ux(end); 

figure(1)
hold on 
h1 = plot(x(1,1:200), (ux_0(1,:)), 'r', 'LineWidth', 2.5)
h2 = plot(x(1,1:200), (uy_0(1,:)), 'b', 'LineWidth', 2.5)

plot(x_fea, ux_fea, 'r--', 'LineWidth', 2.5)
plot(x_fea, uy_fea, 'b--', 'LineWidth', 2.5)
xlabel("X-Coord", 'Interpreter', 'LaTeX')
ylabel("Displacement", 'Interpreter', 'LaTeX')
leg = legend([h1, h2],["X-Disp", "Y-Disp"], 'Location', 'Best');

saveas(gcf, 'BouchbinderDisp.eps', 'epsc')
