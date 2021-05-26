
%% 
alpha = linspace(0,1,100);

a = (1-alpha).^2;
b = (1-alpha).^3;

figure(1)
hold on 
h1 = plot(alpha, a, 'LineWidth', 2.5)
h2 = plot(alpha, b, 'LineWidth', 2.5)
legend([h1, h2], "a(\alpha)", "a^3(\alpha)")
xlabel("Damage")
ylabel("Modulation Functions")

ax = gca; 
ax.FontSize = 25;

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
