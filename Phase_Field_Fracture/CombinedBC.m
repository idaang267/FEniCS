clear all
close all
test_par = 1.0; 
% h0 is the height 
h0 = 10000.0;
% Amount of displacement 
Delta = 10; 
% stretch in direction of displacement 
lambda_a = 1 + Delta/h0;

% Known a for specific geometry 
a = 2*sqrt(h0/pi)*(lambda_a - inv(lambda_a));
c = 1.55;
%% Define LHS section - discontinuous section 
x_l_t = linspace(-0.5, -0.5, 101);
y_l_t = linspace(0.5, 0, 101);
% Convert coordinates
r = sqrt(x_l_t.^2 + y_l_t.^2);
th = atan(y_l_t./x_l_t);

x_l_b = linspace(-0.5, -0.5, 101);
y_l_b = linspace(0, -0.5, 101);
% Convert coordinates
r = sqrt(x_l_b.^2 + y_l_b.^2);
th = atan(y_l_b./x_l_b);
thdeg  = rad2deg(th-pi)

y1_l_t = c*r.*cos(th+pi) - x_l_t;
y2_l_t = a*sqrt(r).*sin((th+pi)./2) - y_l_t;
y1_l_b = c*r.*cos(th-pi) - x_l_b;
y2_l_b = a*sqrt(r).*sin((th-pi)./2) - y_l_b;

%% Define Top section - continuous
x_t_l = linspace(0, -0.5, 101);
y_t_l = linspace(0.5, 0.5, 101);
% Convert coordinates
r = sqrt(x_t_l.^2 + y_t_l.^2);
th = atan(y_t_l./x_t_l);
% th(1) = -th(1);

x_t_r = linspace(0, 0.5, 101);
y_t_r = linspace(0.5, 0.5, 101);
% Convert coordinates
r = sqrt(x_t_r.^2 + y_t_r.^2);
th = atan(y_t_r./x_t_r);

y1_t_l = flip(c*r.*cos(th + pi) - x_t_l);
y2_t_l = flip(a*sqrt(r).*sin((th + pi)./2) - y_t_l);
y1_t_r = c*r.*cos(th) - x_t_r;
y2_t_r = a*sqrt(r).*sin((th)./2) - y_t_r;

%% Define RHS - continuous
x_rhs = linspace(0.5, 0.5, 101);
y_rhs = linspace(0.5, -0.5, 101);
% Convert coordinates
r = sqrt(x_rhs.^2 + y_rhs.^2);
th = atan(y_rhs./x_rhs);
thdeg = rad2deg(th)

y1_rhs = c*r.*cos(th) - x_rhs;
y2_rhs = a*sqrt(r).*sin((th)./2) - y_rhs;

%% Define Bottom - continuous 
x_b_l = linspace(-0.5, 0, 101);
y_b_l = linspace(-0.5, -0.5, 101);
% Convert coordinates
r = sqrt(x_b_l.^2 + y_b_l.^2);
th = atan(y_b_l./x_b_l);
th(101) = -th(101);
thdeg = rad2deg(th-pi)

x_b_r = linspace(0.5, 0, 101);
y_b_r = linspace(-0.5, -0.5, 101);
% Convert coordinates
r = sqrt(x_b_r.^2 + y_b_r.^2);
th = atan(y_b_r./x_b_r);

y1_b_l = flip(c*r.*cos(th-pi) - x_b_l);
y2_b_l = flip(a*sqrt(r).*sin((th-pi)./2) - y_b_l);

y1_b_r = c*r.*cos(th) - x_b_r;
y2_b_r = a*sqrt(r).*sin((th)./2) - y_b_r;

%%
figure(2)
hold on
plot(linspace(-1.5, -1, 101), test_par*y1_b_r, 'r', 'LineWidth', 2.5)
plot(linspace(-1.5, -1, 101), test_par*y2_b_r, 'b', 'LineWidth', 2.5)
plot(linspace(-1, -0.5, 101), test_par*y1_b_l, 'r', 'LineWidth', 2.5)
plot(linspace(-1, -0.5, 101), test_par*y2_b_l, 'b', 'LineWidth', 2.5)
xline(-0.5, '--k')
plot(y_l_t, test_par*y1_l_t, 'r', 'LineWidth', 2.5)
plot(y_l_t, test_par*y2_l_t, 'b', 'LineWidth', 2.5)
h1 = plot(y_l_b, test_par*y1_l_b, 'r', 'LineWidth', 2.5)
h2 = plot(y_l_b, test_par*y2_l_b, 'b', 'LineWidth', 2.5)
xline(0.5, '--k') % top 
plot(linspace(0.5, 1., 101), test_par*y1_t_l, 'r', 'LineWidth', 2.5)
plot(linspace(0.5, 1., 101), test_par*y2_t_l, 'b', 'LineWidth', 2.5)
plot(linspace(1.0, 1.5, 101), test_par*y1_t_r, 'r', 'LineWidth', 2.5)
plot(linspace(1.0, 1.5, 101), test_par*y2_t_r, 'b', 'LineWidth', 2.5)
xline(1.5, '--k') % rhs
plot(linspace(1.5, 2.5, 101), test_par*y1_rhs, 'r', 'LineWidth', 2.5)
plot(linspace(1.5, 2.5, 101), test_par*y2_rhs, 'b', 'LineWidth', 2.5)

leg = legend([h1, h2],["y1", "y2"], 'Location', 'Best');