clear all
close all

% Path to data 
path = "/Users/idaang/Dropbox/cartilage_mechanics_IDA/FEniCS/Swelling/Lagrange/";

% Initiate counter to index each datapoint (indentation and pressure profile)
count = 0;
timestep = [0, 30, 60]; 
%sub = ["1", "1E4", "1E15"];
% SampleRate 
sampleRate = 15;

for m = 1:length(timestep) 
    % Name file consistently 
    %name = [path + num2str(timestep)+ '_' + sub(m) + '.mat'];
    name = [path + num2str(timestep(m)) + '.mat'];
    
    data = load(name);
    Size = size(data.com_Data); 
    line(:,1) = data.com_Data(:,1);
    count = 1; 
    for n = 2:2:Size(2) 
        pressure(:,count) = data.com_Data(:,n); 
        for x = 1:length(pressure)
            if isnan(pressure(x,count))
                indPoint(count,1) = line(x-1,1);
                prePoint(count,1) = -pressure(x-1,count);
                break
            else
                indPoint(count,1) = 1.4;
            end
        end 
        count = count + 1; 
    end 
    
    % Get reference points along X3 
    count2 = 1; 
    for n = 3:2:Size(2)
        refPoint(count2,1) = data.com_Data(1,n);
        count2 = count2 + 1; 
    end
%     figure(1)
%     plot(refPoint,prePoint,'-','LineWidth',2)
%     hold on 
%     xlabel("X_3")
%     ylabel("Pressure")
%     title("Pressure Profile: R = 0.25, d = 0.05 increments")
%     set(gca,'fontsize',12)
    figure(2)
    plot(refPoint,indPoint,'-','LineWidth',2)
    hold on 
    xlabel("X_3")    
    ylabel("u_{top}")
    title("Indentation Profile: R = 0.25, d = 0.05 increments")
    set(gca,'fontsize',12)
end

legend(string(timestep(1)), string(timestep(2)), string(timestep(3)) )
