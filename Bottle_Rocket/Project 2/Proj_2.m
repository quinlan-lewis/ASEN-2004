clear all, close all
%% Constants
g = 9.81; %m/s^2
gamma = 1.4;
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
tspan = [0 5];
g_vec = [0 0 9.8];

%% Initial conditions
P0 = P_gage + P_amb; %pa-inital gage pressure inside bottle
V0 = 0.001; %m^3-initial volume of water inside bottle
T0 = 300; %K-initial air temp
v0 = 0.0; %m/s-initail speed of rocket
th = 45; %initial angle of rocket
x0 = 0.0; %m-initial horixontal distance
y0 = 0.25; %m-initial vertical height
ls = 0.5; %m-length of test stand
M_0 = 1.15;
a_x0 = 0;
a_z0 = 0;
vel_x0 = 0;
vel_z0 = 0;
pos_x0 = 0;
pos_z0 = .25;
Thrusti = 0;


%% Phase 1: Water expulsion

%integrating given equation of force to get velocity

[t,R] = ode45(@(t,R) phase1ODE(t,Cd,At,d_water,P0,gamma,P_amb,Vb,R,M_btl,th,ls,a_x0,a_z0,vel_x0,vel_z0,Ab,g,d_amb,Rey,T0), tspan, [V0,M_0,vel_x0,vel_z0,pos_x0,pos_z0,Thrusti]); %volume of air is V
%tp1 = t
vol_air = zeros(1,65);
Mr = zeros(1,65);
vel_x = zeros(1,65);
vel_z = zeros(1,65);
pos_x = zeros(1,65);
pos_z = zeros(1,65);
Thrust = zeros(1,65);
P_forthrust = zeros(1,65);
    
for i = 1:65
    vol_air(i) = R(i,1);
    Mr(i) = R(i,2);
    vel_x(i) = abs(R(i,3));
    vel_z(i) = abs(R(i,4));
    pos_x(i) = R(i,5);
    pos_z(i) = R(i,6);
end

for i = 1:length(vol_air)
    P_forthrust(i) = ((.001/vol_air(i))^gamma)*P0;
end
for i = 1:length(P_forthrust)
   Thrust(i) = 2*Cd*At*(P_forthrust(i)-P_amb); 
end

%plot of postion for first phase
figure(1)
plot(pos_x,pos_z)
xlim([0,80]);
ylim([0,30]);

%plot of thrust for first phase
figure(2)
plot(Thrust);
ylim([0,200]);
%xlim([0 .45]);


%% Phase 2: Air expulsion
P_end = P0*(V0/Vb)^gamma;
T_end = T0*(V0/Vb)^(gamma-1);

vel_x0 = vel_x(end);
vel_z0 = vel_z(end);
pos_x0 = pos_x(end);
pos_z0 = pos_z(end);

Mr_0 = Mr(end);
Mair_i = (P0*.001)/(Rey*T0);
Thrusti = Thrust(end);
Ve = sqrt(vel_x0^2 + vel_z0^2);

[t,R2] = ode45(@(t,R2) phase2ODE(t,Cd,At,d_water,P0,gamma,P_amb,Vb,R2,M_btl,P_end,T_end,ls,a_x0,a_z0,vel_x0,vel_z0,Ab,g,Rey,Ve,Mr_0,Mair_i), tspan, [vel_x0,vel_z0,pos_x0,pos_z0,Thrusti]);

% for i = 1:65
%     vel_x2(i) = R2(i,1);
%     vel_z2(i) = R2(i,2);
%     pos_x2(i) = R2(i,3);
%     pos_z2(i) = R2(i,4);
% end
hold on figure(1)
plot(pos_x,pos_z)

%Equation 1: sumForces = F_vec-D_vec+M_btl*g_vec
%Equation 22: F = m_air*Ve+(P_amb-P_exit)*At
%Equation 24: M_tot_flow2 = -Cd*d_exit*At*ve



%% Phase 3: Ballisitic 

% vel_x0 = vel_x2(end)
% vel_z0 = vel_z2(end)
% pos_x0 = pos_x2(end)
% pos_z0 = pos_z2(end)
% 
% [t,R3] = ode45(@(t,R3) phase3ODE(t,Cd,At,d_water,P0,gamma,P_amb,Vb,R3,M_btl,P_end,T_end,ls,a_x0,a_z0,vel_x0,vel_z0,Ab,g,Rey,Ve,Mr_0,Mair_i,d_amb), tspan, [vel_x0,vel_z0,pos_x0,pos_z0]);

%Equation 1: sumForces = F_vec-D_vec+M_btl*g_vec
%Equation 25; F = 0, M_tot = M_btl