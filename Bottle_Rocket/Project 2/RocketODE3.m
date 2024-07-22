function Data = RocketODE3(~,Param_0)
global g Cd d_amb Vb P_amb gamma d_water D_thr D_btl At Ab R M_btl CD Ls tspan volAir_0 T0 mAir_0 th P0 Thrusty count

PosX = Param_0(1);
PosZ = Param_0(2);
Vx = Param_0(3);
Vz = Param_0(4);
volAir = Param_0(5);
mAir = Param_0(6);
mRocket = Param_0(7);

V = sqrt(Vx^2 + Vz^2);

D = 0.5*d_amb*(V^2)*CD*Ab;

P_end = P0*(volAir_0/Vb)^gamma;
T_end = T0*(volAir_0/Vb)^(gamma-1);
% fprintf('%f\n',mAir_0)
P_Thrust = P_end*((mAir/mAir_0)^gamma);

d_air_thrust = mAir/Vb;
T_air_thrust = P_Thrust/(d_air_thrust*R);

if Vb > volAir
%     fprintf('P0 = %f\n',P0)
%     fprintf('volAir_0 = %f\n',volAir_0)
%     fprintf('volAir = %f\n',volAir)
%     fprintf('gamma = %f\n',gamma)
    P = P0*((volAir_0/volAir)^gamma);
    
    Ve = sqrt(2*(P-P_amb)/d_water);
    
%     fprintf('Cd = %f\n',Cd)
%     fprintf('At = %f\n',At)
%     fprintf('P = %f\n',P)
%     fprintf('P_amb = %f\n',P_amb)
    F = 2*Cd*At*(P-P_amb);
    Thrusty(count) = F;
    
    dvolAirdt = Cd*At*sqrt((2/d_water)*(P0*((volAir_0/volAir)^gamma)-P_amb));
    
    dmRocketdt = -Cd*At*sqrt(2*d_water*(P-P_amb));
    %dmRocketdt = -Cd*d_water*At*Ve;
    
    dmAirdt = 0;
    
elseif P_Thrust > P_amb
    %P_star = P_Thrust*(2/(gamma+1))^(gamma/(gamma-1));
    P = ((mAir/mAir_0)^gamma)*P_end;
    P_star = P*(2/(gamma+1))^(gamma/(gamma-1));
    T_air_thrust = P/(d_air_thrust*R);
    if P_star > P_amb
        Te = (2/(gamma+1))*T_air_thrust;
        
        Pe = P_star;
        
        de = Pe/(R*Te);
        
        Ve_air = sqrt(gamma*R*Te);
        
        AirFlow = Cd*de*At*Ve_air;
        
        F = AirFlow*Ve_air+(P_amb-Pe)*At;
        Thrusty(count) = F;
    else
        Pe = P_amb;
        
        Me = sqrt((2/(gamma-1))*((P/P_amb)^((gamma-1)/gamma)-1));
        
        Te = T0/(1+((gamma-1)/2)*Me^2);
        
        de = P_amb/(R*Te);
        
        Ve_air = Me*sqrt(gamma*R*Te);
        
        AirFlow = Cd*de*At*Ve_air;
        
        F = AirFlow*Ve_air+(P_amb-Pe)*At;
        Thrusty(count) = F;
    end
    %dmAirdt = -Cd*de*At*Ve_Air;
    dvolAirdt = 0;  
    %F = dmAirdt*VeAir+(P_amb-Pe)*At;
    dmAirdt = -AirFlow;
    dmRocketdt = -AirFlow;
    
else
    F = 0;
    Thrusty(count) = F;
    dmRocketdt = 0;
    dmAirdt = 0;
    dvolAirdt = 0;
end

dxdt = Vx;
dzdt = Vz;

if PosZ < Ls*sind(th)
    Dx = cosd(th);
    Dz = sind(th);
else
    Dx = Vx/V;
    Dz = Vz/V;
end
% fprintf('F = %f\n',F)
% fprintf('D = %f\n',D)
% fprintf('Dx = %f\n',Dx)
% fprintf('mRocket = %f\n',mRocket)
%weight = (mRocket)*9.8;%(d_amb*dvolAirdt)
dVxdt = (((F-D)*Dx)/mRocket);
dVzdt = (((F-D)*Dz)/mRocket)-g;

Data(1) = dxdt;
Data(2) = dzdt;
Data(3) = dVxdt;
Data(4) = dVzdt;
Data(5) = dvolAirdt;
Data(6) = dmAirdt;
Data(7) = dmRocketdt;
    
Data = Data';

count = count+1;
end