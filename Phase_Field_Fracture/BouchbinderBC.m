clear all
close all

%%
% h0 is the height 
h0 = 1.0;
% Amount of displacement 
Delta = 0.01; 
% stretch in direction of displacement 
lambda_a = 1 + Delta/h0;

% Known a for specific geometry 
a = 2*sqrt(h0/pi)*(lambda_a - inv(lambda_a));
c = 2.15;

d_l = 0.005; 

%% Define LHS section - discontinuous section 
x_l_b = linspace(-d_l, -d_l, 101);
y_l_b = linspace(-d_l, -0, 101);
% Convert coordinates
r_l_b = sqrt(x_l_b.^2 + y_l_b.^2);
th_l_b = atan2(y_l_b, x_l_b); % - pi;

x_l_t = linspace(-d_l, -d_l, 101);
y_l_t = linspace(0, d_l, 101);
% Convert coordinates
r_l_t = sqrt(x_l_t.^2 + y_l_t.^2);
th_l_t = atan2(y_l_t, x_l_t); %+ pi;

y1_l_b = c*r_l_b.*cos(th_l_b) - x_l_b;
y2_l_b = a*sqrt(r_l_b).*sin(th_l_b./2) - y_l_b;
y1_l_t = c*r_l_t.*cos(th_l_t) - x_l_t;
y2_l_t = a*sqrt(r_l_t).*sin(th_l_t./2) - y_l_t;

%% Define Top section - continuous
x_t_l = linspace(-d_l, 0, 101);
y_t_l = linspace(d_l, d_l, 101);
% Convert coordinates
r_t_l = sqrt(x_t_l.^2 + y_t_l.^2);
th_t_l = atan2(y_t_l, x_t_l);
% th_t_l(101) = -th_t_l(101);
% th_t_l = th_t_l + pi; 

x_t_r = linspace(0, d_l, 101);
y_t_r = linspace(d_l, d_l, 101);
% Convert coordinates
r_t_r = sqrt(x_t_r.^2 + y_t_r.^2);
th_t_r = atan2(y_t_r, x_t_r);

y1_t_l = c*r_t_l.*cos(th_t_l) - x_t_l;
y2_t_l = a*sqrt(r_t_l).*sin(th_t_l./2) - y_t_l;
y1_t_r = c*r_t_r.*cos(th_t_r) - x_t_r;
y2_t_r = a*sqrt(r_t_r).*sin(th_t_r./2) - y_t_r;

%% Define RHS - continuous
x_rhs = linspace(d_l, d_l, 101);
y_rhs = linspace(d_l, -d_l, 101);
% Convert coordinates
r_r = sqrt(x_rhs.^2 + y_rhs.^2);
th_r = atan(y_rhs./x_rhs);
thdeg = rad2deg(th_r);

y1_rhs = c*r_r.*cos(th_r) - x_rhs;
         %test_par*sqrt(r_r/(2*pi)).*(7/3*cos(th_r./2) - cos(3*th_r./2)) + ...
         %test_par^2*(1/(sqrt(2*pi)))^2*(-1/15*log(r_r) - 52/45*(log(r_r) + 3/4*sin(th_r).^2) - 103/48*cos(th_r) + 26/15*cos(2*th_r) -3/16*cos(3*th_r));
y2_rhs = a*sqrt(r_r).*sin(th_r./2) - y_rhs;
         %test_par*sqrt(r_r/(2*pi)).*(13/3*sin(th_r./2) - sin(3*th_r./2)) + ...
         %test_par^2*(1/(sqrt(2*pi)))^2*(th_r/15 - 52/45*(th_r/4 - 3/8*sin(2*th_r)) - 61/48*sin(th_r) + 26/15*sin(2*th_r) - 3/16*sin(3*th_r)) ;

%% Define Bottom - continuous 
x_b_r = linspace(d_l, 0, 101);
y_b_r = linspace(-d_l, -d_l, 101);
% Convert coordinates
r_b_r = sqrt(x_b_r.^2 + y_b_r.^2);
th_b_r = atan(y_b_r./x_b_r);

x_b_l = linspace(0, -d_l, 101);
y_b_l = linspace(-d_l, -d_l, 101);
% Convert coordinates
r_b_l = sqrt(x_b_l.^2 + y_b_l.^2);
th_b_l = atan(y_b_l./x_b_l);
th_b_l(1) = -th_b_l(1);
th_b_l = th_b_l - pi;

y1_b_r = c*r_b_r.*cos(th_b_r) - x_b_r;
y2_b_r = a*sqrt(r_b_r).*sin(th_b_r./2) - y_b_r;
y1_b_l = c*r_b_l.*cos(th_b_l) - x_b_l;
y2_b_l = a*sqrt(r_b_l).*sin(th_b_l./2) - y_b_l;

%%
figure(2)
hold on
plot(linspace(-d_l*3, -d_l*2, 101), y1_b_r, 'r', 'LineWidth', 2.5)
plot(linspace(-d_l*3, -d_l*2, 101), y2_b_r, 'b', 'LineWidth', 2.5)
plot(linspace(-d_l*2, -d_l, 101), y1_b_l, 'r', 'LineWidth', 2.5)
plot(linspace(-d_l*2, -d_l, 101), y2_b_l, 'b', 'LineWidth', 2.5)
xline(-d_l, '--k')
plot(y_l_t, y1_l_t, 'r', 'LineWidth', 2.5)
plot(y_l_t, y2_l_t, 'b', 'LineWidth', 2.5)
h1 = plot(y_l_b, y1_l_b, 'r', 'LineWidth', 2.5);
h2 = plot(y_l_b, y2_l_b, 'b', 'LineWidth', 2.5);
xline(d_l, '--k') % top 
plot(linspace(d_l, 2*d_l, 101), y1_t_l, 'r', 'LineWidth', 2.5)
plot(linspace(d_l, 2*d_l, 101), y2_t_l, 'b', 'LineWidth', 2.5)
plot(linspace(2*d_l, 3*d_l, 101), y1_t_r, 'r', 'LineWidth', 2.5)
plot(linspace(2*d_l, 3*d_l, 101), y2_t_r, 'b', 'LineWidth', 2.5)
xline(3*d_l, '--k') % rhs
plot(linspace(3*d_l, 5*d_l, 101), y1_rhs, 'r', 'LineWidth', 2.5)
plot(linspace(3*d_l, 5*d_l, 101), y2_rhs, 'b', 'LineWidth', 2.5)

leg = legend([h1, h2],["u1", "u2"], 'Location', 'Best');

%%
r = 1;
th = linspace(-pi, pi, 101);
y1 = c*r.*cos(th);
%   sqrt(r/(2*pi)).*(7/3*cos(th./2) - cos(3*th./2)) + ...
%     (1/(4*sqrt(2*pi)))^2*(-1/15*log(r)-52/45*(log(r)+3/4*sin(th).^2) - 103/48*cos(th) + 26/15*cos(2*th)-3/16*cos(3*th));
y2 = a*sqrt(r).*sin(th./2);

% sqrt(r/(2*pi)).*(13/3*sin(th./2) - sin(3*th./2)) + ...
%     (1/(4*sqrt(2*pi)))^2*(th/15 - 52/45*(th/4 - 3/8*sin(2*th)) - 61/48*sin(th) + 26/15*sin(2*th) - 3/16*sin(3*th)) ;

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