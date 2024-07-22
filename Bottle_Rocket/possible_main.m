%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  %%
%%                   Static Test Stand Data Analysis                      %
%            (I say this title can be edited as we go along - El)         %
% Evan Shults  
% Elena Salgado
%
% Overview:                                                               %
% 1. Determine the specific Impulse (Isp) of the rocket using integration %
%    from teststand data.                                                 %
% 2. Using Isp, determine the deltaV and use as initial condition for the %
%    equations of motion.                                                 %
% 3. Using initial velocity we also assume that the water has been        %
%    expelled, so the mass of the rocket is unchanging and                %
%    is that of the empty bottle. [NOT SURE what this means - Ev]         %
%                                                                         %
% Assuumptions:                                                           %
% - Assumes Instantaneous initial acceleration/velocity                   %
% - High initial drag due to larger velocity                              %
% - Less inertia during initial stages of flight                          %
% - Greater initial angle of trajectory due to large initial velocity     %
% - Simpler incorportaion of statisitcs (mean Isp for model)[NOT SURE-Ev] %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Housekeeping
close all;clear;clc;

%% Constants (Please help, these determine our ISP values)!!!!!!!!!!!!!!!!!!!
%This is to calculate Isp, but idk how to account for the changing mass
%while the water is being expelled, so idk what m should be
%also, idk what the units are for the test stand data, idk what g0 would be
%either

%These are my best guessses and that seems to work out:
m = 0.15; %kg????? possibly 1.1568(from 2012 project)
g0 = 9.81;  %m/s^2????

%% Read in Data Files

files = dir('LA_Demo_30*');


data = {};
for i = 1:length(files)
    data(i) = {readmatrix(files(i).name)};
end


%% %%%%%%%% Convert data from a cell to a table %%%%%%%%%%%%%%%%%

data = cell2table(data);

% Group 301
% water = 1000; %[g]
% pressure = 40; %[psi]
% Temp = 9; %[C]

test1 = data{1,1}{:,:};
test2 = data{1,2}{:,:};
test3 = data{1,3}{:,:};
test4 = data{1,4}{:,:};

% Group 302
% water = 1000; %[g]
% pressure = 40; %[psi]
% Temp = 15; %[C]

test5 = data{1,5}{:,:};
test6 = data{1,6}{:,:};
test7 = data{1,7}{:,:};
test8 = data{1,8}{:,:};
test9 = data{1,9}{:,:};

% Group 303
% water = 1000; %[g]
% pressure = 40; %[psi]
% Temp = %ranges from 9 to 25; %[C]

test10 = data{1,10}{:,:};
test11 = data{1,11}{:,:};
test12 = data{1,12}{:,:};
test13 = data{1,13}{:,:};
test14 = data{1,14}{:,:};
test15 = data{1,15}{:,:};
test16 = data{1,16}{:,:};
test17 = data{1,17}{:,:};
test18 = data{1,18}{:,:};

%put it all in an easier to acces matrix:
test = {test1,test2,test3,test4,test5,test6,test7,test8,test9,test10,test11,test12,test13,test14,test15,test16,test17,test18}';


%% %%%%%%% Remove excess data %%%%%%%%%%%%%
n = 18; %cuase there are 18 graphs


timeidx = zeros(2,18);

for i = 1:n
    
    testi = test{i};
    SummedLoad = testi(:,3);
    
    t = (0:1:(length(SummedLoad)-1))./1000;
    
    
    %Find when the y component is greater than 10 on both sidse of the
    %graph and then subtract and add 0.2 s on the x axis form that index

    idx = (find(SummedLoad > 10));
    start = idx(1)- 100;
    
    %idx = (find(SummedLoad > 10));
    finish = idx(end) + 300;
    
    testi = testi(start:finish,1:3);
    test{i} = testi;
    
    %make sure the time has the same number of elements
    timeidx(:,i) = [start,finish]; 
    
end


%% Mad Plots yo
%in this section we get: plots, Isp, peak thrust, and thrust time
Isp = zeros(1,18);
ThrustTimes = zeros(2,18);
PeakThrust = zeros(1,18);


for i = 1:n  %PLOTS WOO
    %% Sko buffs (-El)
    
    t = timeidx(:,i);       %get the specific time vector 
    t = (t(1):t(2))./1000;
    
    testi = test{i};        %get the specific graph
    name = ['Static Test Integration Curve, Test ',int2str(i)];
    figure()
    hold on
    grid on

    %Plot all three force data lines:
    for j = 1:3
        plot(t,testi(:,j));
    end
    
    
    %% Make xlines:
    %the first xline is when the slopt of the curve is > 1
    %the second xline line when the slope of the curve has toned down
    
    %These xlines represent when the force is applied thus what we should
    %integrae over to find Isp
    
    %first xline:
    SummedLoad = testi(:,3);
    idx = zeros(1,length(SummedLoad)-1);
    for h = 2:length(SummedLoad)
        if (abs(SummedLoad(h)-SummedLoad(h-1)) > 1 && SummedLoad(h) > 0)
            idx(h) = h;
        else
            idx(h) = 0;
        end   
    end
    idx2 = find(idx);
    
    xline(t(idx(idx2(1))),'r');
    thrustStartidx = idx(idx2(1)); %this is to find the duration of thrust
    thrustStart = t(idx(idx2(1)));
    
    %second xline:
    for h = length(SummedLoad):-1:length(SummedLoad)/2
        if abs(SummedLoad(h)-SummedLoad(h-1)) > 0.4
            idx(h) = h;
        else
            idx(h) = 0;
        end   
    end
    idx2 = find(idx);
    
    xline(t(idx(idx2(end))),'r');
    thrustFinidx = idx(idx2(end)); %this is to find the duration of thrust
    thrustFin = t(idx(idx2(end)));
    
    %make yline at 0:
    yline(0,'r');
    
    %other plot stuff:
    title(name)
    xlabel('Time (ms)')
    ylabel('Force (lbf)')
    yline(0);
    ylim([-5 60]);
    legend('Load Cell 1', 'Load Cell 2','Summed Load')
    hold off
    
    %% Find peak thrust and time of thrust
    ThrustTimes(:,i) = [thrustStart;thrustFin]; %this is to find the duration of thrust
    PeakThrust(i) = max(testi(:,j));
    
    %% Now to integrate between those lines (Find Isp)
    intT = (ThrustTimes(1,i):0.001:ThrustTimes(2,i)); %the time to integrate over
    Isp(i) = trapz(intT,testi((thrustStartidx:thrustFinidx),j))/(m*g0); 
    
end

ThrustTime = ThrustTimes(2,:)-ThrustTimes(1,:);
AvgThrustTime = mean(ThrustTime) %Average Thrust Time (s)
AvgPeakThrust = mean(PeakThrust) %Average Peak Thrust (lbf)
AvgIsp = mean(Isp) %Average Isp (s)

%% Make histogrmas of Isp, peak thrust, and total time of thrust for each data set
%%OR we could just take the tables of these values instead if these
%%histograms, I don't have a preference

figure
histogram(Isp,8);
xlabel('Isp Values (s)')
ylabel('# of Occurances')
title('Histogram of Isp Values for Each Data Set')

figure
histogram(PeakThrust,8,'FaceColor','b');
xlabel('Peak Thrust (lbf)')
ylabel('# of Occurances')
title('Histogram of Peak Thrust for Each Data Set')

figure
histogram(ThrustTime,8,'FaceColor','g');
xlabel('Thrust Time (s)')
ylabel('# of Occurances')
title('Histogram of Time of Thrust for Each Data Set')













