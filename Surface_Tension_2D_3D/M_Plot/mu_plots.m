clear all

directory = ["alt_chi_0.4_lm_1.0_mu_1.0_l0_2.50/"];

for i = 1:length(directory)
    dist = importfile(directory(i) + "data_plot.txt",1,inf);
    dataplot = importfile2(directory(i) + "dict.txt",1,inf);
    dataMat = table2array(dataplot); 
    time = dataMat(:,1)';
    d = dist.e2';
    
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
