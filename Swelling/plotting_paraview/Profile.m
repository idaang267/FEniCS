% Path to data 
path = "/Users/idaang/Dropbox/cartilage_mechanics_IDA/FEniCS/Swelling/Data/";

% Initiate counter to index each datapoint (indentation and pressure profile)
count = 0;
% Name samples according to the time step/sample 
timeSample = [0, 10, 20, 29];
% For each time sample 
for m = timeSample
    count = count + 1;
    for i = 1:15
        % Name file consistently 
        name = [path+num2str(m)+'_'+num2str(i)+'.csv'];
        % Read in full csv file
        data = readtable(name);
        % arc_length gives the height of the cube
        len = data.arc_length;
        % Contact Pressure 
        pressure = data.ContactPressure;
        % Current location along the X_3 direction 
        X(count,i) = data.Points_2(1);
        % For every point along the height of the cube
        for j = 1:length(len)
            % Not a Number (NaN) are produced where the cube height is
            % indented.
            if isnan(pressure(j)) 
                % Use NaNs to find the indentation profile 
                pro(count,i) = len(j-1);
                % Find the pressure along the indentation 
                pre(count,i) = -pressure(j-1); 
                % Break when both points are found
                break
            end
        end
    end     
end

%% Indentation Profile Plot
figure(1) 
count = 1; 
for m = timeSample 
    plot(X(count,:),pro(count,:),'-','LineWidth',2)
    hold on 
    count = count + 1; 
end
xlabel("X_3")
ylabel("u_{top}")
title("Indentation Profile: R = 0.25, d = 0.01 increments")
set(gca,'fontsize',12)
legend(string(timeSample(1)), string(timeSample(2)), string(timeSample(3)), string(timeSample(end)))
hold off

%% Pressure Profile Plot
figure(2) 
count = 1;
for m = timeSample
    plot(X(count,:),pre(count,:),'-','LineWidth',2);
    hold on 
    count = count + 1; 
end
xlabel("X_3")
ylabel("Pressure")
title("Pressure Profile: R = 0.25, d = 0.01 increments")
set(gca,'fontsize',12)
legend(string(timeSample(1)), string(timeSample(2)), string(timeSample(3)), string(timeSample(end)))
