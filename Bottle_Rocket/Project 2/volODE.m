function [dvdt] = volODE(t,Cd,At,d_water,P0,gamma,P_amb,Vb,R,M_btl,th,ls,ax_0,az_0,vel_x0,vel_z0,Ab,g)
vol_air = R(1);
if vol_air < Vb
    %% volume of air and mass of rocket
    vol_air = R(1);
    vol_air_dot = Cd*At*sqrt((2/d_water)*(P0*((.001/vol_air)^gamma)-P_amb));
    
    P = ((.001/R(1))^gamma)*P0;
    Ve = sqrt(2*(P-P_amb)/d_water);
    mdot = -Cd*At*(2*d_water*(P-P_amb))^1/2;
    if R(2) > 0.15
        R(2) = mdot;
    else
        R(2) = 0;
    end
    R(1) = vol_air_dot;
    
    %% Thrust
%     R(3) = ax_0;
%     R(4) = az_0;
%     R(5) = vel_x0;
%     R(6) = vel_z0;
    veltot = sqrt(R(5)^2 + R(6)^2);
    Thrust = 2*Cd*At*(P-P_amb);
    
    if ls*cos(th) > .01
        xdir = cos(th);
        zdir = sin(th);
    else 
        xdir = R(5)/veltot;
        zdir = R(6)/veltot;
    end
    
    xThrust = Thrust*sin(atan(zdir/xdir));
    zThrust = Thrust*cos(atan(zdir/xdir));
    
    Drag = -(d_water/2)*(veltot^2)*Cd*Ab;
    xDrag = -(d_water/2)*(R(5)^2)/((cos(atan(zdir/xdir)))^2)*Cd*Ab;
    zDrag = -(d_water/2)*(R(6)^2)/((cos(atan(zdir/xdir)))^2)*Cd*Ab;
    
    accx = ((Thrust+Drag)*xdir)/mdot;
    accz = ((Thrust+Drag)*zdir)/mdot - g;
    R(3) = accx; 
    R(4) = accz; 
    
    velx = Ve*cos(atan(zdir/xdir));
    velz = Ve*sin(atan(zdir/xdir));
    
    R(5) = velx;
    R(6) = velz;
    R(7) = Thrust;
    
else
    R(1) = 0;
    R(2) = 0;
    R(3) = 0;
    R(4) = 0;
    R(5) = 0;
    R(6) = 0;
    R(7) = 0;
end
    
dvdt = R;
    
end

