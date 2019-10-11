%% Access later after saving Mat Files
clear all
close all

% Path to data 
path = "/Users/idaang/Dropbox/cartilage_mechanics_IDA/FEniCS/Swelling/PenaltyGap/";

% Initiate counter to index each datapoint (indentation and pressure profile)
count = 0;
timestep = [1 2]; 
%sub = ["1", "1E4", "1E15"];
% SampleRate 
sampleRate = 30;

for m = 1:length(timestep) 
    % Name file consistently 
    %name = [path + num2str(timestep)+ '_' + sub(m) + '.mat'];
    name = [path + num2str(timestep(m)) + '.mat'];
    
    data = load(name);
    Size = size(data.com_Data); 
    line(:,1) = data.com_Data(:,1);
    count = 1; 
    for n = 3:3:Size(2) 
        Gap(:,count) = data.com_Data(:,n); 
        for x = 1:length(Gap)
            if isnan(Gap(x,count))
                indPoint(count,1) = line(x-1,1);
                prePoint(count,1) = Gap(x-1,count);
                break
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
    
    figure(1)
    plot(refPoint,prePoint,'-','LineWidth',2)
    hold on 
    xlabel("X_3")
    ylabel("Gap")
    title("R = 0.25, d = 0.05 increments")
    set(gca,'fontsize',20)
    figure(2)
    plot(refPoint,indPoint,'-','LineWidth',3)
    hold on 
    xlabel("X_3")    
    ylabel("u_{top}")
    title("R = 0.25, d = 0.05 increments")
    set(gca,'fontsize',20)
end

legend(string(timestep(1)), string(timestep(2)), string(timestep(3)), string(timestep(4)))
%legend(string(timestep(1)), string(timestep(2)), string(timestep(3)), string(timestep(4)), string(timestep(5)), string(timestep(6)), string(timestep(end)))