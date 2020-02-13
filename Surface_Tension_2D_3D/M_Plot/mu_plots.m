clear all

directory = ["alt_chi_0.2_l0_2.00/"];

for i = 1:length(directory)
    dist = importfile(directory(i) + "data_plot.txt",1,inf);
    dataplot = importfile2(directory(i) + "dict.txt",1,inf);
    dataMat = table2array(dataplot); 
    time = dataMat(:,1);
    
    [x, y] = size(dataMat); 
    chemPotCoord = [];
 
    for a = 1:x
        count1 = 1;
        count2 = 1;
        for b = 2:y
            if mod(b,2) == 0 
                Coord(count1,a) = dataMat(a,b);  
                count1 = count1 + 1;
            else 
                chemPotCoord(count2,a) = dataMat(a,b);
                count2 = count2 + 1; 
            end 
        end        
    end 
end

%%
n = 10^-3;
gamma = 1.0;
chi = 0.4;
l0 = 3.0;
J0 = l0^3;
mu0 = n*(1/l0-1/J0) + 1/J0 + log(1 - 1/J0) + chi/(J0^2)+ n*gamma*2/l0

