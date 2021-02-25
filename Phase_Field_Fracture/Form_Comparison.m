
%% 
alpha = linspace(0,1,100);

a = (1-alpha).^2;
b = (1-alpha).^3;

figure(2)
plot(alpha, b, "--")
alpha_t = 0.7;

for i = 1:length(b)
    if alpha_t <= alpha(i)
        b(i) = 0;
    end   
end

figure(2)
hold on 
title("Criterion: alpha_t = 0.7")
h1 = plot(alpha, a, 'LineWidth', 2.5)
h2 = plot(alpha, b, 'LineWidth', 2.5)
legend([h1, h2], "a(alpha)", "b(alpha)")
xlabel("Damage")

%% New 
% undamaged (0) to damaged (1)
alpha = linspace(0,1,100);

a = (1-alpha).^2;
b = (1-alpha).^3;

figure(3)
plot(alpha, b, "--")

alpha_t = 0.1;

for i = 1:length(b)
    if alpha_t >= alpha(i)
        b(i) = 0;
    end   
end

figure(3)
hold on 
title("Criterion switched: alpha_t = 0.1")
h1 = plot(alpha, a, 'LineWidth', 2.5)
h2 = plot(alpha, b, 'LineWidth', 2.5)
legend([h1, h2], "a(alpha)", "b(alpha)")
xlabel("Damage")
