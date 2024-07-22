function [dvdt] = phase3ODE(t,Cd,At,d_water,P0,gamma,P_amb,Vb,R3,M_btl,P_end,T_end,ls,a_x0,a_z0,vel_x0,vel_z0,Ab,g,Rey,Ve,Mr_0,Mair_i,d_amb)
    velx = R3(3)
    velz = R3(4)
    veltot = sqrt((velx)^2+(velz)^2);
    Drag = -0.5 * d_amb * Cd * Ab * veltot^2;
    mass = Mr_0;

    dirx = velx/veltot;
    dirz = velz/veltot;

    Thrust = 0;

    accelx = ((Thrust + Drag)*dirx)/mass;
    accelz = ((Thrust + Drag)*dirz)/mass - g;

    dvdt = [accelx,accelz,velx,velz]';
        
end