%% Housekeeping
close all
clear
clc

%%%%%%%%%%%%%%% DESCRIPTION %%%%%%%%%%%%%%%%
%ikeMain reads in the raw data files and runs the ISPcalc function to
%return the ISP (ISP) max thrust (Tpeak), and time elapsed (t)

%Date Created: 4/3/20
%Date Modified: 4/12/20

%Author: Isaak Timko
%SID: 107 668 395

%%%%%%%%%%% DATA{ 1:4 } %%%%%%%%%%%%
% Group 301
% water     = 1000;     %[g]
% pressure  = 40;       %[psi]
% Temp      = 9;        %[C]

%%%%%%%%%%% DATA{ 5:9 } %%%%%%%%%%%%
% Group 302
% water     = 1000;     %[g]
% pressure  = 40;       %[psi]
% Temp      = 15;       %[C]

%%%%%%%%%%% DATA{ 10:18 } %%%%%%%%%%%%
% Group 303
% water     = 1000;     %[g]
% pressure  = 40;       %[psi]
% Temp      = 9:25;     %[C]

%% Read in Data Files

files   = dir('LA_Demo_30*');
n = length(files);

data    = cell(1,n);
for i = 1:n
    data{1,i}   = readmatrix(files(i).name);
end

data{1}     = data{1}(3000:3700,1:3); 
data{2}     = data{2}(2550:3250,1:3);
data{3}     = data{3}(2250:2950,1:3);
data{4}     = data{4}(3000:3700,1:3);
data{5}     = data{5}(2850:3550,1:3);
data{6}     = data{6}(2750:3450,1:3); 
data{7}     = data{7}(850:1550,1:3);
data{8}     = data{8}(1550:2250,1:3);
data{9}     = data{9}(2600:3300,1:3);
data{10}    = data{10}(2800:3500,1:3);
data{11}    = data{11}(2250:2950,1:3);
data{12}    = data{12}(2950:3650,1:3);
data{13}    = data{13}(2950:3650,1:3);
data{14}    = data{14}(2300:3000,1:3);
data{15}    = data{15}(4500:5200,1:3);
data{16}    = data{16}(3150:3850,1:3);
data{17}    = data{17}(3400:4100,1:3);
data{18}    = data{18}(2700:3400,1:3);

%% Converting from lbf to N

for i = 1:n
    data{i} = data{i} *4.44822; %N
end

%% Calculating Isp

ISP     = zeros(1,n);
Tpeak   = zeros(1,n);
t       = zeros(1,n);

for i = 1:n
    [ISP(i), Tpeak(i), t(i)] = ISPcalc(data{i});
end

%% Plotting Test stand data

%plotting the force for each individual test
% for i = 1:n
%    figure
%    x = linspace(0,t(i),length(data{i}(:,3)));
%    hold on
%    plot(x,data{i}(:,1))
%    plot(x,data{i}(:,2))
%    plot(x,data{i}(:,3))
%    hold off
%    title(strcat('Test #',num2str(i)))
%    xlabel('Time (s)')
%    ylabel('Thrust (N)')
%    legend('Transducer 1','Transducer 2','Total Force','Location','Best')
% end

%plotting the total force for all 18 tests
figure
hold on
elements = cell(1,n);
for i = 1:n
   x = linspace(0,t(i),length(data{i}(:,3)));
   plot(x,data{i}(:,3))
   elements{i} = strcat('Test #', num2str(i));
end
title('All Tests')
xlabel('Time (s)')
ylabel('Thrust (N)')
legend(elements,'Location','Best')

%% Error Analysis

% Calculating The mean values
MeanISP     = mean(ISP);    %seconds
MeanTpeak   = mean(Tpeak);  %N
MeanTime    = mean(t);      %seconds

%calculating the Standard Deviation
StdISP      = std(ISP);     %seconds
StdTpeak    = std(Tpeak);   %N
StdTime     = std(t);       %seconds

%calculating the standard error of the mean
SEMISP      = StdISP/sqrt(length(ISP));
SEMTpeak    = StdTpeak/sqrt(length(Tpeak));
SEMTime     = StdTime/sqrt(length(t));

%Confidence Interval 
z           = [1.96 2.24 2.58]; % 95% 97.5% 99%
CIISP       = z*SEMISP;
CITpeak     = z*SEMTpeak;
CITime      = z*SEMTpeak;

%Computing SEM with increasing tests
SEMVec          = zeros(1,n);
for i = 1:n
    StdTEMP     = std(ISP(1:i));
    SEMVec(i)   = StdTEMP/sqrt(i);
end
SEMVec(1) = SEMVec(2) + (SEMVec(2)-SEMVec(3));

%Computing SEM projected for many tests

reps    = 100;
SEMProj = zeros(1,reps-n);

for i = n:reps
    SEMProj((i-n)+1)  = StdISP/sqrt(i);
end

% figure
% hold on
% plot(SEMVec)
% plot(n:reps,SEMProj)
% hold off
% title('SEM vs. N')
% xlabel('Number Of Trials')
% ylabel('SEM')
% legend('Provided Data','Projected Data','Location','best')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% New code????
global g Cd d_amb Vb P_amb gamma d_water D_thr D_btl At Ab R M_btl CD Ls tspan volAir_0 T0 mAir_0 th P0 vq

interpdata = abs(data{18}(:,3) * 4.448);
time = linspace(0,t(1),length(data{18}(:,3)));
time_ellapsed = t(length(t));
xq = 0:(time(2)-time(1))/4:time(length(time));
tspan = [0 time(length(time))];

%Constants
g = 9.81; %m/s^2
gamma = 1.4;
Cd = 0.8; %discharge coeficient
d_amb = 0.961; %kg/m^3-ambient air density
Vb = 0.002; %m^3-volume of empty bottle
P_amb = 12.1*6894.75729; %pa-atmospheric pressure
P_gage = 72.5*6894.75729; %pa- gage pressure
d_water = 1000; %kg/m^3-density of water
D_thr = 2.1/100; %m-diameter of throat
D_btl = 10.5/100; %m-diameter of bottle
R = 287; %J/Kg*K-gas constant of air
M_btl = 0.15; %Mass of empty rocket
CD = 0.5; %drag coeficient
At = pi*(D_thr/2)^2; %Throat Area
Ab = pi*(D_btl/2)^2; %bottle area
Ls = 0.5; %m-length of test stand

%% Initial conditions
P0 = P_gage + P_amb; %pa-inital gage pressure inside bottle
V0 = 0.001; %m^3-initial volume of water inside bottle
T0 = 300; %K-initial air temp
v0 = 0.0; %m/s-initail speed of rocket
th = 45; %initial angle of rocket
x0 = 0.0; %m-initial horixontal distance
z0 = 0.25; %m-initial vertical height

%assigning initial condition values
M_H20 = d_water*V0;
v_x0 = 0;
v_z0 = 0;
volAir_0 = 0.001;
mAir_0 = (P0*V0)/(R*T0);
mRocket_0 = M_btl + M_H20 + mAir_0;

vq = interp1(time,interpdata,xq);
figure()
hold on
plot(time,interpdata,'o',xq,vq,':.');

[t,mass] = ode45('InterpODE',tspan,vq);
    
% accel = vq/mass;
% [t,vel] = ode45('velODE',tspan,accel);
% [t,pos] = ode45('posODE',tspan,vel);
% 
% figure()
% plot(pos)




