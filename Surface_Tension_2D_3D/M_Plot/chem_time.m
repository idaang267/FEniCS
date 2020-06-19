clear all
close all

format shortEng

directory = ["exp_chi_0.2_g_1.00_l0_2.00_exp_1.6/"];

load('Pore_Case3_ET_p_dict.mat')

Chi = 0.2; 
NOmega = 10^-3; 
Lambda0 = 2.00; 
gamma = 1.00; 
J0 = Lambda0^3;

%Time = [1E-6, 1E-3, 1E-2, 0.1, 1.0];

mu0 = NOmega*(1/Lambda0-1/J0) + 1/J0 + log(1 - 1/J0) + Chi/(J0^2) + (NOmega*gamma*2)/Lambda0;
NuEq = 0.5 - NOmega/2*(1/(Lambda0^2*(Lambda0^3-1)) + NOmega/Lambda0^2 - 2*Chi/Lambda0^5 )^(-1);
TimeNorm = Time_list.*(1/NOmega)*(Lambda0^6/(Lambda0^3 - 1))*(3*(1+2*NuEq)/(2*(3+5*NuEq)));

% Black, Blue, Green, Purple, light blue, red, dark red
colors = {[0 0 0], [0 0.4470 0.7410], [0.4 0.8 0.2],...
          [0.4940 0.1840 0.5560], [0.3010 0.7450 0.9330],...
          [0.9 0.01 0.01], [0.6350 0.0780 0.1840]};

colRes = 1; 
for i = 1:length(directory)   
    dist = importfile2(directory(i) + "rad_chem_time.txt",1,inf);
    dataMat = table2array(dist); 
    time = dataMat(1,:);
    finRad = dataMat(2,:);
    rad = dataMat(3:end,1); 
    count = 1;
    
    for a = 2+10:length(time)
        if isnan(time(a)) == true
            lastReadTime = a;
            break
        end
    end 
    
    for a = 2:lastReadTime-1
        chem_pot = dataMat(3:end,a);
        
        figure(1)
        grid on 
        set(gca,'fontsize',15)
        xlabel("X")
        ylabel("$\tilde{\mu}$", 'Interpreter', 'Latex')
        plot(rad,chem_pot,'-','Color', colors{colRes},'LineWidth',2)
        hold on   
        
        curRad(:,a) = rad.*finRad(a); 
        maxCurRadArr(a) = curRad(end,end);
        
        figure(2) 
        grid on
        set(gca,'fontsize',15)
        xlabel("x")
        ylabel("$\tilde{\mu}$", 'Interpreter', 'Latex')
        plot(curRad(:,a),chem_pot,'-','Color', colors{colRes},'LineWidth',2)
        hold on
        maxCurrRad = max(maxCurRadArr);
        
        % For legend
        newTime{count} = num2str(time(a),'%10.2e\n');
        count = count + 1;
        if colRes == length(colors)
            colRes = 1; 
        else
            colRes = colRes + 1; 
        end
    end
end 

figure(1)
leg = legend(newTime, 'Location', 'Best');  
title(leg,'Time');
xlim([0 1])
%ylim([-5.5E-3 0])

figure(2)
leg = legend(newTime, 'Location', 'Best');  
title(leg,'Time');
xlim([0,maxCurrRad])
%ylim([-5.5E-3 0])

%% 
for i = 1:length(Time_list) 
   legTime{i} = num2str(TimeNorm(i), '%10.2e\n');
end 

Time0Norm = - tID0*NOmega/Lambda0 + mu0;
Time1Norm = - tID1*NOmega/Lambda0 + mu0;
Time2Norm = - tID2*NOmega/Lambda0 + mu0;
Time3Norm = - tID3*NOmega/Lambda0 + mu0;
Time4Norm = - tID4*NOmega/Lambda0 + mu0; 
Time5Norm = - tID5*NOmega/Lambda0 + mu0; 

figure(5)
hold on 
plot(R,Time0Norm,'-.','LineWidth',2)
plot(R,Time1Norm,'-.','LineWidth',2)
plot(R,Time2Norm,'-.','LineWidth',2)
plot(R,Time3Norm,'-.','LineWidth',2)
plot(R,Time4Norm,'-.','LineWidth',2)
plot(R,Time5Norm,'-.','LineWidth',2)
hold off
xlabel("X")
ylabel("$\tilde{\mu}$", 'Interpreter', 'Latex')
leg = legend({legTime{1}, legTime{2}, legTime{3}, legTime{4}, legTime{5}, legTime{6}}, 'Location', 'Best', 'NumColumns',2);
title(leg,'Time');
xlim([0 1])
set(gca,'fontsize', 22)

%%
format long
NewTime(1) = time(31)*NOmega*((Lambda0^3 - 1)/Lambda0^6)*(2*(3+5*NuEq)/(3*(1+2*NuEq)))
NewTime(2) = time(35)*NOmega*((Lambda0^3 - 1)/Lambda0^6)*(2*(3+5*NuEq)/(3*(1+2*NuEq)))
NewTime(3) = time(38)*NOmega*((Lambda0^3 - 1)/Lambda0^6)*(2*(3+5*NuEq)/(3*(1+2*NuEq)))
NewTime(4) = time(40)*NOmega*((Lambda0^3 - 1)/Lambda0^6)*(2*(3+5*NuEq)/(3*(1+2*NuEq)))
NewTime(5) = time(42)*NOmega*((Lambda0^3 - 1)/Lambda0^6)*(2*(3+5*NuEq)/(3*(1+2*NuEq)))
NewTime(6) = time(46)*NOmega*((Lambda0^3 - 1)/Lambda0^6)*(2*(3+5*NuEq)/(3*(1+2*NuEq)))

CheckTimeNorm = NewTime.*(1/NOmega)*(Lambda0^6/(Lambda0^3 - 1))*(3*(1+2*NuEq)/(2*(3+5*NuEq)));

%%

%saveas(gcf, 'RadCur_Simp_Chem_chi_0.2_g_1.00_l0_2.00.eps', 'epsc')
saveas(gcf, "swelling_radius_norm_option2.fig")
% saveas(gcf, 'Legend2.eps', 'epsc')
% saveas(gcf, "Option2.fig")
