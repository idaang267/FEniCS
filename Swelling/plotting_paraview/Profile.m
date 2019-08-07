path = "/Users/idaang/Dropbox/cartilage_mechanics_IDA/FEniCS/Swelling/Data/";

count = 0;
timeSample = [10,20,25,45];
for m = timeSample
    count = count + 1;
    for i = 1:15
        name = [path+num2str(m)+'_'+num2str(i)+'.csv'];
        data = readtable(name);
        Chem = data.ChemicalPotential;
        len = data.arc_length;
    
        X(count,i) = data.Points_2(1);
        for j = 1:length(len)
            if isnan(Chem(j)) 
                pro(count,i) = len(j);
                break
            end
        end
    end
    plot(X(count,:),pro(count,:),'-','LineWidth',2)
    hold on 
end

xlabel("X_3")
ylabel("u_{top}")
title("R = 0.175, d = 0.01 increments")
set(gca,'fontsize',12)

legend(string(timeSample(1)), string(timeSample(2)), string(timeSample(3)),string(timeSample(end)))
