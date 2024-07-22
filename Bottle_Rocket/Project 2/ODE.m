function calc = ODE(~,Initials)
    global g gamma Cd d_amb Vb P_amb P_gage d_water D_thr D_btl Rey M_btl CD At Ab tspan P0 T0
    
    x0 = Initials(1);
    z0 = Initials(2);
    vx0 = Initials(3);
    vz0 = Initials(4);
    vol_air0 = Initials(5);
    m_air0 = Initials(6);
    m_rocket0 = Initials(7);

    veltot = sqrt(vx0^2+vz0^2);
    
    P_end = P0*(vol_air0/Vb)^gamma;
    T_end = T0*(vol_air0/Vb)^(gamma-1);
    
    
    
end

