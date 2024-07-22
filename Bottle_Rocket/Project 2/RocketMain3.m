clear all, close all
%% Constants
global g Cd d_amb Vb P_amb gamma d_water D_thr D_btl At Ab R M_btl CD Ls tspan volAir_0 T0 mAir_0 th P0 Thrusty count

g = 9.81; %m/s^2
gamma = 1.4;
Cd = 0.8; %discharge coeficient
d_amb = 0.961; %kg/m^3-ambient air density
Vb = 0.002; %m^3-volume of empty bottle
P_amb = 12.1*6894.75729; %pa-atmospheric pressure
P_gage = 50*6894.75729; %pa- gage pressure
d_water = 1000; %kg/m^3-density of water
D_thr = 2.1/100; %m-diameter of throat
D_btl = 10.5/100; %m-diameter of bottle
R = 287; %J/Kg*K-gas constant of air
M_btl = 0.15; %Mass of empty rocket
CD = 0.5; %drag coeficient
At = pi*(D_thr/2)^2; %Throat Area
Ab = pi*(D_btl/2)^2; %bottle area
tspan = [0 5];
Ls = 0.5; %m-length of test stand
Thrusty = [1,113];
count = 1;
%% Initial conditions
P0 = P_gage + P_amb; %pa-inital gage pressure inside bottle
V0 = 0.001; %m^3-initial volume of water inside bottle
T0 = 300; %K-initial air temp
v0 = 0.0; %m/s-initail speed of rocket
th = 45; %initial angle of rocket
x0 = 0.0; %m-initial horixontal distance
z0 = 0.25; %m-initial vertical height

M_H20 = d_water*V0;
v_x0 = 0;
v_z0 = 0;
volAir_0 = 0.001;
mAir_0 = (P0*V0)/(R*T0);
mRocket_0 = M_btl + M_H20 + mAir_0;

Param_0 = [x0,z0,v_x0,v_z0,volAir_0,mAir_0,mRocket_0];
%[t,data] = ode45(@(t,Data) RocketODE2(t,Param_0), tspan, Param_0);
[t,data] = ode45('RocketODE3', tspan, Param_0);
Xpos = data(:,1);
Zpos = data(:,2);
Vx = data(:,3);
Vz = data(:,4);
volAir = data(:,5);
mAir = data(:,6);
mRocket = data(:,7);

figure(1)
plot(Xpos,Zpos)
xlim([0,100]);
ylim([0,30]);
Pair = zeros(1,length(Xpos));
Thrust = zeros(1,length(Xpos));
time = linspace(0,0.5,length(Xpos));
for i = 1:length(Xpos)
    if volAir(i) < Vb
    PairW(i) = P0*(volAir_0/volAir(i))^gamma;
    Thrust(i) = 2*Cd*At*(PairW(i)-P_amb);
    else
        Pair(i) = PairW(end)*(mAir(i)/mAir_0)^gamma;
        if Pair(i) > P_amb
            P_star(i) = Pair(i)*(2/(gamma+1))^(gamma/(gamma-1));
            rhoe(i) = mAir(i)/Vb;
            T(i) = Pair(i)/(rhoe(i)*R);
            if P_star > P_amb
                               
                Te(i) = (2/(gamma+1))*T(i);
                
               
                Pe(i) = P_star(i);
                rhoe(i) = Pe(i)/(R*Te(i));
                
                Ve(i) = sqrt(gamma*R*Te(i));
                
                Thrust(i) = Cd*rhoe(i)*At*((Ve(i))^2)+(P_amb-Pe(i))*At;
            else
                Pe(i) = P_amb;
                
                Me(i) = sqrt((2/(gamma-1))*((Pair(i)/P_amb)^((gamma-1)/gamma)-1));
                
                Te(i) = T0/(1+((gamma-1)/2)*Me(i)^2);
                
                Ve(i) = Me(i)*sqrt(gamma*R*Te(i));
                
                Thrust(i) = Cd*rhoe(i)*At*((Ve(i))^2)+(P_amb-Pe(i))*At;
            end
        end
                
    end
    
        
end

plot(time,Thrust)
