Points = DJDataS3.Points0;
J002 = DJDataS3.Y002; 
J008 = DJDataS3.Y008;
D002 = DJDataS3.Y1;
D008 = DJDataS3.Y2; 

%%
hold on 
yyaxis left 
h1 = plot(Points, J002, 'LineWidth', 2.5)
plot(Points, ones(length(Points)), 'k:')
ylabel('J')
yyaxis right
h3 = plot(Points, D002, 'LineWidth', 2.5) 
ylabel('Damage')
xlabel('X-Coordinates')
ax = gca; 
ax.FontSize = 20;
xlim([0, 0.25])


%%
hold on 
yyaxis left 
h1 = plot(Points, J002, 'LineWidth', 2.5)
h2 = plot(Points, J008, 'LineWidth', 2.5)
% plot(Points, ones(length(Points)), 'k:')
ylabel('J')
yyaxis right
h3 = plot(Points, D002, 'LineWidth', 2.5) 
h4 = plot(Points, D008, 'LineWidth', 2.5) 
legend([h1, h3, h2, h4], "Y = 0.02", "Y = 0.02", "Y = 0.05", "Y = 0.05")
ylabel('Damage')
xlabel('X-Coordinates')
ax = gca; 
ax.FontSize = 20;
xlim([0, 0.25])

%%
hold on 
h1 = plot(Points, J002, 'LineWidth', 2.5)
h2 = plot(Points, J008, 'LineWidth', 2.5)
plot(Points, ones(length(Points)), 'k:')
ylabel('J')
legend([h1, h2], "Y = 0.02", "Y = 0.05")
xlabel('X-Coordinates')
ax = gca; 
ax.FontSize = 20;
xlim([0, 0.25])

