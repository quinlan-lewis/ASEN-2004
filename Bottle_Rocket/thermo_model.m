clear,close all

%Author: Quinlan Lewis
%Date: 4/20/20
%Student ID's: 108502480

%**THIS CODE WAS MODIFIED FROM ASEN 2012 PROJECT 2 CODE**%

%% Constants
%initializing all of the given values for the different constants
%creating globals to be used in our ode function
global g Cd d_amb Vb P_amb gamma d_water D_thr D_btl At Ab R M_btl CD Ls tspan volAir_0 T0 mAir_0 th P0 v_windy v_windx v_windz

g = 9.81; %m/s^2
gamma = 1.4;
Cd = 0.8; %discharge coeficient
d_amb = 0.961; %kg/m^3-ambient air density
Vb = 0.002; %m^3-volume of empty bottle
P_amb = 14.45156*6894.75729; %pa-atmospheric pressure
P_gage = 40.5*6894.75729; %pa- gage pressure
d_water = 1000; %kg/m^3-density of water
D_thr = 2.1/100; %m-diameter of throat
D_btl = 10.5/100; %m-diameter of bottle
R = 287; %J/Kg*K-gas constant of air
M_btl = 121/1000; %Mass of empty rocket
CD = 0.327; %drag coeficient
At = pi*(D_thr/2)^2; %Throat Area
Ab = pi*(D_btl/2)^2; %bottle area
tspan = [0 5]; %time of bottle rocket flight
Ls = 0.5; %m-length of test stand

%% Initial conditions
P0 = P_gage + P_amb; %pa-inital gage pressure inside bottle
V0 = 0.001; %m^3-initial volume of water inside bottle
T0 = 22+273; %K-initial air temp
v0 = 0.0; %m/s-initail speed of rocket
th = 45; %initial angle of rocket
x0 = 0.0; %m-initial horixontal distance
z0 = 0.25; %m-initial vertical height
y0 = 0;

%assigning initial condition values
% M_H20 = d_water*V0;
M_H20 = 993/1000;
v_windy = (1*.44704); %m/s wind velocity
v_windz = 0;
v_windx = 0;
v_x0 = 0;
v_z0 = 0;
v_y0 = 0;
volAir_0 = 0.001;
mAir_0 = (P0*V0)/(R*T0);
mRocket_0 = M_btl + M_H20 + mAir_0;

%creating array of initial values
Param_0 = [x0,z0,y0,v_x0,v_z0,v_y0,volAir_0,mAir_0,mRocket_0];

opts=odeset('Events',@stopping_point);
[t, X] = ode45(@thermo_modelODE, [0 10], Param_0, opts);
tspan = [0 t(end)];

%call ode45 function to numerically integrate given equations and result in
%values for position in the x and z, velocity in x and z, volume of air in
%bottle, mass of air in the bottle, and mass of the rocket over the flight
%of the rocket
[time,data] = ode45('thermo_modelODE', tspan, Param_0);%returns the data of from the ode45 function
fatThrust = zeros(1,length(time));%initializes the length of the thrust for the rocket

%below we assign the values returned from the RocketODE call based off of
%the entered initial conditions
Xpos = data(:,1);
Zpos = data(:,2);
Ypos = data(:,3);
Vx = data(:,4);
Vz = data(:,5);
Vy = data(:,6);
volAir = data(:,7);
mAir = data(:,8);
mRocket = data(:,9);

%now that all of the data is sorted into appropriate arrays we can index
%them to graph them
figure(1)
plot3(Xpos,Ypos,Zpos)
grid on
xlim([0,80]);
zlim([0,25]);
title('Plot of Gold Rockets trajectory')
xlabel('meters')
ylabel('meters')
zlabel('meters')

x_data = zeros(100,1);
y_data = zeros(100,1);

for i = 1:100
    th = 44 + (46-44).*rand(1,1);
    M_H20 = 992.5/1000 + (993.5/1000-992.5/1000).*rand(1,1);
    mRocket_0 = M_btl + M_H20 + mAir_0;
    v_windy = (0.5*.44704) + ((1.5-0.5)*.44704).*rand(1,1);
    P_gage = (39*6894.75729) + ((41-39)*6894.75729).*rand(1,1);
    P0 = P_gage + P_amb;
    
    Param_0 = [x0,z0,y0,v_x0,v_z0,v_y0,volAir_0,mAir_0,mRocket_0];
    
    opts=odeset('Events',@stopping_point);
    [t, X] = ode45(@thermo_modelODE, [0 10], Param_0, opts);
    tspan = [0 t(end)];
    
    [time,data] = ode45('thermo_modelODE', tspan, Param_0);
    x_data(i) = data(length(data),1);
    y_data(i) = data(length(data),3);
end
 
%% Replace this section of code with your real data
% Simulate and plot 100 data points, (YOU SHOULD USE REAL DATA HERE!)
N = 100; % Number of points to simulate
%%

figure(6); subplot(2,1,1); plot(x_data,y_data,'k.','markersize',6,'DisplayName','Monte-Carlo Sim points')
axis equal; grid on; xlabel('x [m]'); ylabel('y [m]'); hold on;
 
% Calculate covariance matrix
P = cov(x_data,y_data);
mean_x = mean(x_data);
mean_y = mean(y_data);
 
% Calculate the define the error ellipses
n=100; % Number of points around ellipse
p=0:pi/n:2*pi; % angles around a circle
 
[eigvec,eigval] = eig(P); % Compute eigen-stuff
xy_vect = [cos(p'),sin(p')] * sqrt(eigval) * eigvec'; % Transformation
x_vect = xy_vect(:,1);
y_vect = xy_vect(:,2);
 
% Plot the error ellipses overlaid on the same figure
plot(1*x_vect+mean_x, 1*y_vect+mean_y, 'b','DisplayName','1 Standard Deviation')
plot(2*x_vect+mean_x, 2*y_vect+mean_y, 'g','DisplayName','2 Standard Deviations')
plot(3*x_vect+mean_x, 3*y_vect+mean_y, 'r','DisplayName','3 Standard Deviations')
plot(mean_x,mean_y,'c*','markersize',8,'DisplayName','Prediction Point')
legend
title('Error Ellipses and Monte-Carlo sim data')
%str = 'Prediction Point';
%text(mean_x,mean_y,str)

%% Sensitivity Analysis
%restating constants and initial conditions
g = 9.81; %m/s^2
gamma = 1.4;
Cd = 0.8; %discharge coeficient
d_amb = 0.961; %kg/m^3-ambient air density
Vb = 0.002; %m^3-volume of empty bottle
P_amb = 14.45156*6894.75729; %pa-atmospheric pressure
P_gage = 40*6894.75729; %pa- gage pressure
d_water = 1000; %kg/m^3-density of water
D_thr = 2.1/100; %m-diameter of throat
D_btl = 10.5/100; %m-diameter of bottle
R = 287; %J/Kg*K-gas constant of air
M_btl = 126/1000; %Mass of empty rocket
CD = 0.327; %drag coeficient
At = pi*(D_thr/2)^2; %Throat Area
Ab = pi*(D_btl/2)^2; %bottle area
tspan = [0 5]; %time of bottle rocket flight
Ls = 0.5; %m-length of test stand

%% Initial conditions
P0 = P_gage + P_amb; %pa-inital gage pressure inside bottle
V0 = 0.001; %m^3-initial volume of water inside bottle
T0 = 10+273; %K-initial air temp
v0 = 0.0; %m/s-initail speed of rocket
th = 45; %initial angle of rocket
x0 = 0.0; %m-initial horixontal distance
z0 = 0.25; %m-initial vertical height
y0 = 0;

%assigning initial condition values
% M_H20 = d_water*V0;
M_H20 = 983/1000;
v_windy = (7*.44704); %m/s wind velocity
v_windz = 0;
v_windx = 0;
v_x0 = 0;
v_z0 = 0;
v_y0 = 0;
volAir_0 = 0.001;
mAir_0 = (P0*V0)/(R*T0);
mRocket_0 = M_btl + M_H20 + mAir_0;
%Vary the follow parameters while keeping the others constant:
%1. Coefficient of drag
%2. Mass of the water propellant
%3. Density of the water propellant
%4. Temperature of the water propellant
%5. Launch pad angle

%varying coefficient of drag
CD_vec = 0.2:.01:0.4; %drag coeficient
M_H20 = 983/1000;
mRocket_0 = M_btl + M_H20 + mAir_0;
d_water = 1000; %kg/m^3-density of water
th = 45; %initial angle of rocket
figure(); hold on; grid on
for i = 1:length(CD_vec)
    CD = CD_vec(i);
    
    Param_0 = [x0,z0,y0,v_x0,v_z0,v_y0,volAir_0,mAir_0,mRocket_0];
    
    opts=odeset('Events',@stopping_point);
    [t, X] = ode45(@thermo_modelODE, [0 10], Param_0, opts);
    tspan = [0 t(end)];
    
    [time,data] = ode45('thermo_modelODE', tspan, Param_0);
    x_data = data(end,1);
    scatter(CD,x_data,'b');
end
title('Varying The Drag Coefficient')
ylabel('Downrange Distance(m)')
xlabel('Coefficient of Drag')

%varying mass of propelant
CD_vec = 0.327; %drag coeficient
M_H20_vec = (900:1:1100)/1000;
d_water = 1000; %kg/m^3-density of water
th = 45; %initial angle of rocket
figure(); hold on; grid on
for i = 1:length(M_H20_vec)
    M_H20 = M_H20_vec(i);
    mRocket_0 = M_btl + M_H20 + mAir_0;
    
    Param_0 = [x0,z0,y0,v_x0,v_z0,v_y0,volAir_0,mAir_0,mRocket_0];
    
    opts=odeset('Events',@stopping_point);
    [t, X] = ode45(@thermo_modelODE, [0 10], Param_0, opts);
    tspan = [0 t(end)];
    
    [time,data] = ode45('thermo_modelODE', tspan, Param_0);
    x_data = data(end,1);
    scatter(M_H20,x_data,'b')
end
title('Varying Mass of Propellant')
xlabel('Mass of Propellant(kg)')
ylabel('Downrange Distance(m)')

%varying density of water
CD = 0.327; %drag coeficient
M_H20 = 983/1000;
mRocket_0 = M_btl + M_H20 + mAir_0;
d_water_vec = (1:.001:1.06)*1000; %kg/m^3-density of water
th = 45; %initial angle of rocket
figure(); hold on; grid on
for i = 1:length(d_water_vec)
    d_water = d_water_vec(i);
%     M_H20 = d_water*V0;
%     mRocket_0 = M_btl + M_H20 + mAir_0;
    
    Param_0 = [x0,z0,y0,v_x0,v_z0,v_y0,volAir_0,mAir_0,mRocket_0];
    
    opts=odeset('Events',@stopping_point);
    [t, X] = ode45(@thermo_modelODE, [0 10], Param_0, opts);
    tspan = [0 t(end)];
    
    [time,data] = ode45('thermo_modelODE', tspan, Param_0);
    x_data = data(end,1);
    scatter(d_water,x_data,'b')
end
title('Varying Density of Water')
xlabel('Density of Water(kg/m^3)')
ylabel('Downrange Distance(m)')

%varying temperature of water propellant
CD = 0.327; %drag coeficient
M_H20 = 983/1000;
mRocket_0 = M_btl + M_H20 + mAir_0;
T_water = 0:5:100;
d_water_vec = [.99987,1,.99973,.99913,.99823,.99707,.99567,.994,.99224,.990,.98807,.985,.98324,.9811,.97781,.9752,.97183,.9683,.964,.961,.95838]*1000; %kg/m^3-density of water
%^ courtesy of: https://www.researchgate.net/figure/Density-and-Surface-Tension-of-Water-Against-Air-at-Various-Temperatures_tbl1_290544699
th_vec = 45; %initial angle of rocket
figure(); hold on; grid on
for i = 1:length(T_water)
    d_water = d_water_vec(i);
    
    Param_0 = [x0,z0,y0,v_x0,v_z0,v_y0,volAir_0,mAir_0,mRocket_0];
    
    opts=odeset('Events',@stopping_point);
    [t, X] = ode45(@thermo_modelODE, [0 10], Param_0, opts);
    tspan = [0 t(end)];
    
    [time,data] = ode45('thermo_modelODE', tspan, Param_0);
    x_data = data(end,1);
    scatter(T_water(i),x_data,'b')
end
title('Varying Water Temperature')
xlabel('Temperature(deg C)')
ylabel('Downrange Distance(m)')

%varying launch angle of the pad
CD = 0.327; %drag coeficient
M_H20 = 983/1000;
mRocket_0 = M_btl + M_H20 + mAir_0;
d_water = 1000; %kg/m^3-density of water
th_vec = 35:1:90; %initial angle of rocket
figure(); hold on; grid on
for i = 1:length(th_vec)
    th = th_vec(i);
    
    Param_0 = [x0,z0,y0,v_x0,v_z0,v_y0,volAir_0,mAir_0,mRocket_0];
    
    opts=odeset('Events',@stopping_point);
    [t, X] = ode45(@thermo_modelODE, [0 10], Param_0, opts);
    tspan = [0 t(end)];
    
    [time,data] = ode45('thermo_modelODE', tspan, Param_0);
    x_data = data(end,1);
    scatter(th,x_data,'b')
end
title('Varying Launch Angle')
xlabel('Launch Angle(degrees)')
ylabel('Downrange Distance(m)')

%% Error Analysis
%lower Drag Coef. resulted in most increase range
x_data = zeros(100,1);
y_data = zeros(100,1);
CD = 0.21;

for i = 1:100
    T0 = (9.5+273) + ((10.5-9.5)+273).*rand(1,1);
    th = 44.5 + (45.5-44.5).*rand(1,1);
    M_H20 = 982.5/1000 + (983.5/1000-982.5/1000).*rand(1,1);
    M_btl = 125.5/1000 + ((126.5 - 125.5)/1000).*rand(1,1);
    mRocket_0 = M_btl + M_H20 + mAir_0;
    v_windy = (6*.44704) + ((8-6)*.44704).*rand(1,1);
    P_gage = (39*6894.75729) + ((41-39)*6894.75729).*rand(1,1);
    P0 = P_gage + P_amb;
    
    Param_0 = [x0,z0,y0,v_x0,v_z0,v_y0,volAir_0,mAir_0,mRocket_0];
    
    opts=odeset('Events',@stopping_point);
    [t, X] = ode45(@thermo_modelODE, [0 10], Param_0, opts);
    tspan = [0 t(end)];
    
    [time,data] = ode45('thermo_modelODE', tspan, Param_0);
    x_data(i) = data(length(data),1);
    y_data(i) = data(length(data),3);
end 
N = 100; % Number of points to simulate

figure(6); subplot(2,1,2); plot(x_data,y_data,'k.','markersize',6,'DisplayName','Monte-Carlo Sim points')
axis equal; grid on; xlabel('x [m]'); ylabel('y [m]'); hold on;
 
% Calculate covariance matrix
P = cov(x_data,y_data);
mean_x = mean(x_data);
mean_y = mean(y_data);
 
% Calculate the define the error ellipses
n=100; % Number of points around ellipse
p=0:pi/n:2*pi; % angles around a circle
 
[eigvec,eigval] = eig(P); % Compute eigen-stuff
xy_vect = [cos(p'),sin(p')] * sqrt(eigval) * eigvec'; % Transformation
x_vect = xy_vect(:,1);
y_vect = xy_vect(:,2);
 
% Plot the error ellipses overlaid on the same figure
plot(1*x_vect+mean_x, 1*y_vect+mean_y, 'b','DisplayName','1 Standard Deviation')
plot(2*x_vect+mean_x, 2*y_vect+mean_y, 'g','DisplayName','2 Standard Deviations')
plot(3*x_vect+mean_x, 3*y_vect+mean_y, 'r','DisplayName','3 Standard Deviations')
plot(mean_x,mean_y,'c*','markersize',8,'DisplayName','Prediction Point')
legend
title('Error Ellipses and Monte-Carlo sim data, CD = 0.21')
%str = 'Prediction Point';
%text(mean_x,mean_y,str)