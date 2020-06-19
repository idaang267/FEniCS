clear all

directory = ["exp_chi_0.2_g_1.00_l0_2.00_exp_1.6/"];

for i = 1:length(directory)
    dataplot = importfile2(directory(i) + "dict.txt",1,inf);
    dataMat = table2array(dataplot); 
    time = dataMat(:,1)';
    dist = importfile(directory(i) + "data_plot.txt",1,inf);
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
