%Author: August Hauter
%Author SID: 109116865

%Student Names:
%1) August Hauter 
%2) Tianyi Pu

function fxns = diffEQs (~,Initials)
%% Bringing in variables & Initial Housekeeping
global g Cd rhoAir VolBottle Pamb gamma rhoWater DThroat DBottle AThroat ...
    ABottle R MBottle CD Lstand tspan Volair0 Tair0 Mair0 theta Pair0

x = Initials(1);
z = Initials(2);
vx = Initials(3);
vz = Initials(4);
VolAir = Initials(5);
mAir = Initials(6);
Mrocket = Initials(7);

%Drag Equation
%Setting up velocity for drag equation
vMag = sqrt(vx^2+vz^2); %[m/s]

D = 0.5*rhoAir*(vMag^2)*CD*ABottle; %[N]

%Defining the Pressure and Temperature of the air after the all of the
%water is gone
Pend = Pair0*(Volair0/VolBottle)^gamma; %[Pa] P right after all the water is gone
Tend = Tair0*(Volair0/VolBottle)^(gamma-1); %[K] T right after all the water is gone

PairThrust = Pend*(mAir/Mair0)^gamma; %[Pa] pressure of the air in the bottle 
                                      %during the air thrust phase
rhoAirThrust = mAir/VolBottle; %[kg/m^3] density of the air in the bottle during
                               %the air thrust phase
TairThrust = PairThrust/(rhoAirThrust*R); %[K] temperature of the air in the bottle
                                          %during the air thrust phase

%% Setting up equations for each phase of flight
if VolAir < VolBottle % Setting up differential equations for ode45 Before Water is Exhausted
    %Setting up equation for total pressure of air inside the bottle
    Pair = Pair0*(Volair0/VolAir)^gamma; %[Pa]

    %Setting up the mass flow rate of water out of the bottle
    ve = sqrt(2*(Pair-Pamb)/rhoWater); %[m/s] velocity of water as it exits the bottle

    massFlowWater = Cd*rhoWater*AThroat*ve; %[kg/s]

    %Thrust Equation
    F = 2*Cd*AThroat*(Pair-Pamb); %[N] 
    
    %Setting up the change in air volume over time
    dVolAirdt = Cd*AThroat*sqrt((2/rhoWater)*(Pair0*((Volair0/VolAir)^(gamma))-Pamb)); %[m^3/s]
    dMrocketdt = -massFlowWater; %[kg/s] Change in mass over time for whole rocket
    dmAirdt = 0; %[kg/s] mass flow of the air out of the bottle

elseif PairThrust > Pamb  % Setting up differential equations for ode45 After Water is Exhausted
    %Defining the critical Pressure to determine the exit velocity of the air
    Pstar = PairThrust*(2/(gamma+1))^(gamma/(gamma-1)); %[Pa]

    if Pstar >= Pamb %Choked flow M=1
        Te = (2/(gamma+1))*TairThrust; %[K] temperature of air at exit
        Pe = Pstar; %[Pa] pressure of air at exit
        rhoe = Pe/(R*Te); %[kg/m^3] density of air at exit
        veAir = sqrt(gamma*R*Te); %[m/s] velocity of the air exiting the bottle
    
        massFlowAir = Cd*rhoe*AThroat*veAir; %[kg/s] mass flow of the air out of the bottle
        F = massFlowAir*veAir+(Pamb-Pe)*AThroat; %[N] Thrust from the air
    else %Non-choked flow M<1
        Me = sqrt((2*((PairThrust/Pamb)^((gamma-1)/(gamma)))-2)/(gamma-1)); %[] exit Mach of air
        Te = Tair0/(1+(gamma-1)*Me^2/2); %[K] temperature of the air exiting the 
                                        %bottle during the air thrust phase
        Pe = Pamb; %[Pa] Pressure of the air exiting the bottle
        rhoe = Pamb/(R*Te); %[kg/m^3] density of the air exiting the bottle
        veAir = Me*sqrt(gamma*R*Te); %[m/s] velocity of the air exiting the bottle
    
        massFlowAir = Cd*rhoe*AThroat*veAir; %[kg/s] mass flow of the air out of the bottle
        F = massFlowAir*veAir+(Pamb-Pe)*AThroat; %[N] Thrust from the air
    end
    dVolAirdt = 0;
    dmAirdt = -massFlowAir;
    dMrocketdt = -massFlowAir;
else %Setting up the differential equations for ode45 for the Ballistic Phase of flight
    F = 0;
    dVolAirdt = 0;
    dmAirdt = 0;
    dMrocketdt = 0;
end
%% Setting Motion ODEs to be solved

%ODEs of the position over time
dxdt = vx; %[m/s] derivative of the x position vs time
dzdt = vz; %[m/s] derivative of the z position vs time

%Heading vector
if x < Lstand*cosd(theta)
    Hx = cosd(theta);
    Hz = sind(theta);
else
    Hx = dxdt/vMag; %[m/s] unit vector of the x-velocity
    Hz = dzdt/vMag; %[m/s] unit vector of the z-velocity
end

%Acceleration vector based off the sum of the forces acting on the rocket
dvxdt = (F-D)*(Hx/Mrocket); %[N] the horizontal acceleration of the rocket
dvzdt = (F-D)*(Hz/Mrocket)-g; %[N] the vertical acceleration of the rocket

%% Outputing fxns to be inputted into ode45
fxns(1) = dxdt;
fxns(2) = dzdt;
fxns(3) = dvxdt;
fxns(4) = dvzdt;
fxns(5) = dVolAirdt;
fxns(6) = dmAirdt;
fxns(7) = dMrocketdt;
%ode45 needs a column vector therefore:
fxns = fxns'; 
end