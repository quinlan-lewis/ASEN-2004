clear; clc; close all

%% Constants
%initializing all of the given values for the different constants
%creating globals to be used in our ode function
global g Cd d_amb Vb P_amb gamma d_water D_thr D_btl At Ab R M_btl CD Ls tspan volAir_0 T0 mAir_0 th P0

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
tspan = [0 5]; %time of bottle rocket flight
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

%creating array of initial values
Param_0 = [x0,z0,v_x0,v_z0,volAir_0,mAir_0,mRocket_0];

%call ode45 function to numerically integrate given equations and result in
%values for position in the x and z, velocity in x and z, volume of air in
%bottle, mass of air in the bottle, and mass of the rocket over the flight
%of the rocket
[t,data] = ode45('RocketODE2', tspan, Param_0);%returns the data of from the ode45 function
fatThrust = zeros(1,length(t));%initializes the length of the thrust for the rocket

%The below for loop runs through the RocketODE function to find the thrust
%at each iteration of the integration, because our RocketODE function
%returns an array of thrust and differential equations we can have it
%return only the thrust rather than the data from the differential
%equations
for i = 1:length(t)
   [~,fatThrust(i,1)] = RocketODE2(t(i),data(i,:)); 
end

%below we assign the values returned from the RocketODE call based off of
%the entered initial conditions
Xpos = data(:,1);
Zpos = data(:,2);
Vx = data(:,3);
Vz = data(:,4);
volAir = data(:,5);
mAir = data(:,6);
mRocket = data(:,7);

%now that all of the data is sorted into appropriate arrays we can index
%them to graph them
figure(1)
plot(Xpos,Zpos)
xlim([0,100]);
ylim([0,50]);
hold on
%This for loop finds where phase 1 of the rocket ends where all of the water is expelled and returns the index
%of where that happens as the variable x
for i = 1:length(t)
   if volAir(i) > Vb
       x = i;
       break;
   end
end
%using this index we can graph a dot as a place holder at that index
p1 = plot(Xpos(x),Zpos(x),'.');
p1.MarkerSize = 20;

%This for loop finds where phase 2 of the flight ends and the rocket goes
%into the ballistic phase where there is no thrust and returns that index
%as y
for i = 1:length(t)
   if fatThrust(i) == 0
       y = i;
       break;
   end
end
%using this index we can plot the point at which the rocket enters the next
%phase of its flight
p2 = plot(Xpos(y),Zpos(y),'.');
p2.MarkerSize = 20;
xlabel('Distance in x direction [m]');
ylabel('Height of rocket in z direction [m]')
title('Trajectory of Rocket')