%Authors: Quinlan Lewis, Patrick Tippens
%Date: 12/10/19
%Student ID's: 108502480, 109209805

function [Data, fatThrust] = RocketODE2(~,Param_0)
%% Entering globals to be used in the rest of the function and taking in and assigning the initial conditions
global g Cd d_amb Vb P_amb gamma d_water D_thr D_btl At Ab R M_btl CD Ls tspan volAir_0 T0 mAir_0 th P0 

PosX = Param_0(1);%initial x position [m]
PosZ = Param_0(2);%initial z position [m]
Vx = Param_0(3);%initial x velocity [m/s]
Vz = Param_0(4);%initial v velocity [m/s]
volAir = Param_0(5);%initial volume of air [m^3]
mAir = Param_0(6);%initial mass of the air in the rocket [kg]
mRocket = Param_0(7);%initial mass of the rocket [kg]

V = sqrt(Vx^2 + Vz^2);%magnitude of velocity [m/s]

D = 0.5*d_amb*(V^2)*CD*Ab;%drag force [N]

P_end = P0*(volAir_0/Vb)^gamma;%Pressure of bottle after all of water has been expelled [Pa]
T_end = T0*(volAir_0/Vb)^(gamma-1);%Temperature of bottle after all of water has been expelled [Pa]
P_Thrust = P_end*((mAir/mAir_0)^gamma);%Pressure during the second phase of flight [Pa]

d_air_thrust = mAir/Vb;%density of the air during phase 2 [kg/m^3]
T_air_thrust = P_Thrust/(d_air_thrust*R);%temperature of air during phase 2 [K]

%% Phase 1 of Rocket flight
if Vb > volAir %check to see if the rocket is still in phase 1 by looking at the volume of the air
    P = P0*((volAir_0/volAir)^gamma);%calculate the pressure from the current volume of air [Pa]
    
    Ve = sqrt(2*(P-P_amb)/d_water);%exit velocity of water [m/s]
    
    F = 2*Cd*At*(P-P_amb);%thrust of the rocket in phase 1 [N]
    
    %differential equation for change of volume of air in the bottle [m^3]
    VolAirChange = Cd*At*sqrt((2/d_water)*(P0*((volAir_0/volAir)^gamma)-P_amb));
    
    %differential equation for change in mass of rocket [kg]
    MassRocketChange = -Cd*At*sqrt(2*d_water*(P-P_amb));
    
    %differential equation for change of mass of air in the bottle [kg]
    MassAirChange = 0;
    
%% Phase 2 of Rocket flight
elseif P_Thrust > P_amb %check to see if the rocket is in phase 2
    
    P = ((mAir/mAir_0)^gamma)*P_end;%calculate pressure for phase 2 [Pa]
    P_star = P*(2/(gamma+1))^(gamma/(gamma-1));%calculate critical pressure for phase2 [Pa]
    T_air_thrust = P/(d_air_thrust*R);%Temperature of air at exit during phase 2 [K]
    
    if P_star > P_amb %check to see if flow is choked
        
        Te = (2/(gamma+1))*T_air_thrust;%Calculate the exit temperature [K]
        Pe = P_star;%assign exit pressure [Pa]        
        de = Pe/(R*Te);%calculate the exit density of the air [kg/m^3]
        Ve_air = sqrt(gamma*R*Te);%calculate the exit velocity of the air [m/s]
        
    else %runs if flow is not choked
        
        Pe = P_amb;%Exit pressure [Pa]
        Me = sqrt((2/(gamma-1))*((P/P_amb)^((gamma-1)/gamma)-1));%Mach number at exit
        Te = T0/(1+((gamma-1)/2)*Me^2);%Temperature at exit [K]
        de = P_amb/(R*Te);%density at exit [kg/m^3]
        Ve_air = Me*sqrt(gamma*R*Te);%exit velocity of air [m/s]
        
    end
    %assigning the different differntial equations to the correct variables
    MassAirChange = -Cd*de*At*Ve_air;%given equation for change of mass of air
    VolAirChange = 0;%change of volume of air is 0 because volume of air is now constant
    F = -MassAirChange*Ve_air+(Pe-P_amb)*At;%given equation for thrust
    MassRocketChange = MassAirChange;%given equation for change of mass of rocket
else
    %if not in phase 1 or 2 then it is in phase 3, the ballistic phase
    F = 0;%thrust is 0 because there is no more expulsion
    MassRocketChange = 0;%change of mass of the rocket is 0 because mass is now constant
    MassAirChange = 0;%change of mass of air is 0 because it is now constant
    VolAirChange = 0;%change of volume of air is 0 because the volume of air is now constant
end

%assign the current value of thrust to the return array of thrust
fatThrust = F;

%change in x and z is equal to the velocity in the x and z of the rocket
Velx = Vx;
Velz = Vz;

%Determine if the rocket is still on the stand by checking to see if the
%angle is still constant
if PosZ < Ls*sind(th)
    %assign direction of x and z if it is on the stand
    Dx = cosd(th);
    Dz = sind(th);
else
    %assign direction of x and z if it is no longer on the stand 
    Dx = Vx/V;
    Dz = Vz/V;
end

%calculate the acceleration of the rocket in the x ans z direction subtract
%gravity from the z direction because it acts in the negative z direction
Accelx = (((F-D)*Dx)/mRocket);
Accelz = (((F-D)*Dz)/mRocket)-g;

%assign all of the data to the correct indices to be returned to the main
%function and transpose the ending array because the ODE function requires
%consistent array shapes
Data(1) = Velx;
Data(2) = Velz;
Data(3) = Accelx;
Data(4) = Accelz;
Data(5) = VolAirChange;
Data(6) = MassAirChange;
Data(7) = MassRocketChange;
Data = Data';

end