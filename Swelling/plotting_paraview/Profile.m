%%
clear all
close all

% Path to data 
path = "/Users/idaang/Dropbox/cartilage_mechanics_IDA/FEniCS/Swelling/HE_Data/";

% Initiate counter to index each datapoint (indentation and pressure profile)
count = 0;
% Name samples according to the time step/sample 
timeSample = [32];
%sub = ["1E15"];
% SampleRate 
sampleRate = 45;
% Another counter
total = 1;
% For each time sample 
for m = timeSample
    count = count + 1;
    for i = 1:sampleRate
        % Name file consistently 
        %name = [path+num2str(m)+ '_' + num2str(i) + '_' + sub(count) + '.csv'];
        name = [path + num2str(m)+ '_' + num2str(i) + '.csv'];
        % Read in full csv file
        data = readtable(name);
        % arc_length gives the height of the cube
        len = data.arc_length;
        % Contact Pressure 
        pressure = data.ContactPressure;
        % Save data
        if total == 1
           com_Data(:, total) = len;
        end
        total = total + 1; 
        com_Data(:, total) = pressure;
        total = total + 1;
        com_Data(:, total) = data.Points_2;
    end     
end

%% 

% Initiate counter to index each datapoint (indentation and pressure profile)
count = 0;
% Name samples according to the time step/sample 
% SampleRate 
Dim = size(com_Data);
count = 1;
x = 1;
for i = 2:2:Dim(2)
    for j = 1:Dim(1)
        if isnan(com_Data(j,i)) 
            if mod(count,sampleRate+1) == 0
                x = x + 1;
                count = 1;
            end
            com_pre(count,x) = -com_Data(j-1, i);
            com_pro(count,x) = com_Data(j-1,1);
            com_ind(count,x) = com_Data(j-1,i+1);
            count = count + 1;
            break
        end
    end
end

%% Indentation Profile Plot
figure(1) 
for i = 1:length(timeSample)
    plot(com_ind(:,i),com_pro(:,i),'-','LineWidth',2)
    hold on 
end
xlabel("X_3")
ylabel("u_{top}")
title("Indentation Profile: R = 0.25, d = 0.01 increments")
set(gca,'fontsize',12)
legend(string(sub(1)), string(sub(2)), string(sub(end)))
hold off

%% Pressure Profile Plot
figure(2) 
for i = 1:length(timeSample)
    plot(com_ind(:,i),com_pre(:,i),'-','LineWidth',2)
    hold on 
end
xlabel("X_3")
ylabel("Pressure")
title("Pressure Profile: R = 0.25, d = 0.01 increments")
set(gca,'fontsize',12)
legend(string(sub(1)), string(sub(2)), string(sub(end)))

%% Access later after saving Mat Files
clear all
close all

% Path to data 
path = "/Users/idaang/Dropbox/cartilage_mechanics_IDA/FEniCS/Swelling/HE_Data/";

% Initiate counter to index each datapoint (indentation and pressure profile)
count = 0;
timestep = [0, 10, 20, 30, 32]; 
%sub = ["1", "1E4", "1E15"];
% SampleRate 
sampleRate = 45;

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
    ylabel("Pressure")
    title("Pressure Profile: R = 0.25, d = 0.01 increments")
    set(gca,'fontsize',12)
    figure(2)
    plot(refPoint,indPoint,'-','LineWidth',2)
    hold on 
    xlabel("X_3")    
    ylabel("u_{top}")
    title("Indentation Profile: R = 0.25, d = 0.01 increments")
    set(gca,'fontsize',12)
end

legend(string(timestep(1)), string(timestep(2)), string(timestep(3)), string(timestep(4)), string(timestep(end)))