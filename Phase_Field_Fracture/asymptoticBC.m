%% Correct Cartesian for a square domain 
clear all
close all

r = 1;
th = linspace(-pi, pi, 101);
y1 = r*cos(th);
y2 = sqrt(r)*sin(th./2);

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

%% Define LHS section - discontinuous section 
x_top = linspace(-0.5, -0.5, 101);
y_top = linspace(0.5, 0, 101);
% Convert coordinates
r = sqrt(x_top.^2 + y_top.^2);
th_test = rad2deg(atan(y_top./x_top));
% Theta 
th = linspace(3*pi/4, pi, 101);

y1_top = r.*cos(th);
y2_top = sqrt(r).*sin((th)./2);

figure(2)
hold on
plot(y_top, y1_top, 'r', 'LineWidth', 2.5)
plot(y_top, y2_top, 'b','LineWidth', 2.5)

x_top = linspace(-0.5, -0.5, 101);
y_top = linspace(0, -0.5, 101);
% Convert coordinates
r = sqrt(x_top.^2 + y_top.^2);
th = linspace(-pi, -3*pi/4, 101);

y1_top = r.*cos(th);
y2_top = sqrt(r).*sin((th)./2);

figure(2)
hold on
h1 = plot(y_top, y1_top, 'r', 'LineWidth', 2.5)
h2 = plot(y_top, y2_top,'b' , 'LineWidth', 2.5)
leg = legend([h1, h2],["y1", "y2"], 'Location', 'Best');

%%
x_left_top = linspace(-0.5, -0.5, 101);
y_left = linspace(0.5, 0, 101); 
th_f = rad2deg(atan(y_left./x_left_top) + pi);
th = linspace(3*pi/4, pi, 101);
figure
hold on
plot(y_rhs, th_f, 'LineWidth', 2.5)
plot(y_rhs, rad2deg(th), 'LineWidth', 2.5)

% Test bottom left
x_left_bot = linspace(-0.5, -0.5, 101);
y_left = linspace(-0.5, 0, 101); 
th_f = rad2deg(atan(y_left./x_left_bot));
th = linspace(-pi, -3*pi/4, 101);
figure
hold on
plot(y_rhs, th_f, 'LineWidth', 2.5)
plot(y_rhs, rad2deg(th), 'LineWidth', 2.5)

%% Define Top 
x_top = linspace(0.5, -0.5, 101);
y_top = linspace(0.5, 0.5, 101);
th_f = rad2deg(atan(y_top./x_top));

% Convert coordinates
r = sqrt(x_top.^2 + y_top.^2);
% Theta 
th = linspace(pi/4, 3*pi/4, 101);

y1_top = r.*cos(th);
y2_top = sqrt(r).*sin((th+pi)./2);

figure(3)
hold on
plot(x_top, y1_top, 'r', 'LineWidth', 2.5)
plot(x_top, y2_top, 'b','LineWidth', 2.5)

% figure
% plot(x_top, th_f, 'LineWidth', 2.5)

%% Test top right 
x_top_right = linspace(0.5, 0, 101);
y_top = linspace(0.5, 0.5, 101); 
th_f = rad2deg(atan(y_top./x_top_right));
th = linspace(pi/4, pi/2, 101);
figure
hold on
plot(x_top_right, th_f, 'LineWidth', 2.5)
plot(x_top_right, rad2deg(th), 'LineWidth', 2.5)

% Test top left
x_top_left = linspace(0, -0.5, 101);
y_top = linspace(0.5, 0.5, 101); 
th_f = rad2deg(atan(y_top./x_top_left) + pi);
th = linspace(pi/2, 3*pi/4, 101);
figure
hold on
plot(x_top_left, th_f, 'LineWidth', 2.5)
plot(x_top_left, rad2deg(th), 'LineWidth', 2.5)

%% Define Bottom 
a = 0.014890;
c = 0.011616;
x_bot = linspace(-0.5, 0.5, 101);
y_bot = linspace(-0.5, -0.5, 101);
% Convert coordinates
r = sqrt(x_bot.^2 + y_bot.^2);
th = linspace(-3*pi/4, -pi/4, 101);

y1_bot = r.*cos(th);
y2_bot = sqrt(r).*sin(th./2);

figure(4)
hold on 
plot(x_bot, y1_bot, 'LineWidth', 2.5) 
plot(x_bot, y2_bot, 'LineWidth', 2.5)
th_f = atan(y_bot./(x_bot));

%%
% Test bottom right 
x_bot_right = linspace(0, 0.5, 101);
y_bot = linspace(-0.5, -0.5, 101); 
% y_rhs = linspace(-0.5, 0.5, 101);
th_f = rad2deg(atan(y_bot./x_bot_right));
th = linspace(-pi/2, -pi/4, 101);
figure
hold on
plot(y_rhs, th_f, 'LineWidth', 2.5)
plot(y_rhs, rad2deg(th), 'LineWidth', 2.5)

% Test bottom left
x_bot_right = linspace(-0.5, 0, 101);
y_bot = linspace(-0.5, -0.5, 101); 
% y_rhs = linspace(-0.5, 0.5, 101);
th_f = rad2deg(atan(y_bot./x_bot_right) - pi);
th = linspace(-3*pi/4, -pi/2, 101);
figure
hold on
plot(y_rhs, th_f, 'LineWidth', 2.5)
plot(y_rhs, rad2deg(th), 'LineWidth', 2.5)


%% Define RHS 
x_rhs = linspace(0.5, 0.5, 101);
y_rhs = linspace(-0.5, 0.5, 101);
% Convert coordinates
r = sqrt(x_rhs.^2 + y_rhs.^2);
th = linspace(-pi/4, pi/4, 101);

a = 0.148899;
c = 0.116162;
y1_rhs = c*r.*cos(th);
y2_rhs = a*sqrt(r).*sin(th./2);

figure(5)
hold on
plot(y_rhs, y1_rhs, 'LineWidth', 2.5) 
plot(y_rhs, y2_rhs, 'LineWidth', 2.5)

%%
% Test section
x_rhs = linspace(0.5, 0.5, 101);
y_rhs = linspace(-0.5, 0.5, 101);
th_f = rad2deg(atan(y_rhs./x_rhs));
figure
hold on
plot(y_rhs, th_f, 'LineWidth', 2.5)
plot(y_rhs, rad2deg(th), 'LineWidth', 2.5)
