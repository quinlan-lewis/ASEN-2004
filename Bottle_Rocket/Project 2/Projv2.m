clear all, close all
%% Constants
global g gamma Cd d_amb Vb P_amb P_gage d_water D_thr D_btl Rey M_btl CD At Ab tspan P0 T0
g = 9.81; %m/s^2
gamma = 1.4; %gamma used for air vel calculations
Cd = 0.8; %discharge coeficient
d_amb = 0.961; %kg/m^3-ambient air density
Vb = 0.002; %m^3-volume of empty bottle
P_amb = 12.1*6894.75729; %pa-atmospheric pressure
P_gage = 50*6894.75729; %pa- gage pressure
d_water = 1000; %kg/m^3-density of water
D_thr = 0.021; %m-diameter of throat
D_btl = 0.105; %m-diameter of bottle
Rey = 287; %J/Kg*K-gas constant of air
M_btl = 0.15; %Mass of empty rocket
CD = 0.5; %drag coeficient
At = pi*(D_thr/2)^2; %Throat Area
Ab = pi*(D_btl/2)^2; %bottle area
tspan = [0 5]; %time of rocket launch overall

%% Initial conditions
P0 = P_gage + P_amb; %pa-inital gage pressure inside bottle
T0 = 300; %K-initial air temp
vtot0 = 0.0; %m/s-initail speed of rocket
vx0 = 0;
vz0 = 0; 
th = 45; %initial angle of rocket
x0 = 0.0; %m-initial horixontal distance
y0 = 0.25; %m-initial vertical height
z0 = .01; %m-intial z height
ls = 0.5; %m-length of test stand
vol_air0 = .001;
vol_water0 = .001;
m_water0 = vol_water0*d_water;
m_air0 = (P0*vol_air0)/(Rey*T0);
m_rocket0 = m_air0 + m_water0 + M_btl;

Intials = [x0 z0 vx0 vz0 vol_air0 m_air0 m_rocket0];

[t,Results] = ode45(@(t,Results) ODE(t,Initials),tspan,Initials);