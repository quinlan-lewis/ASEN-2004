function [dmdt] = mODE(t,Cd,At,d_water,P0,gamma,P_amb,M_btl)
mrocket = M_btl+1;
if mrocket > M_btl
    P = ((.001/vol_air)^gamma)*P0;
    Ve = sqrt(2*(P-P_amb)/d_water);
    %mdot = -Cd*d_water*At*Ve;
    mdot = -Cd*At*(2*d_water*(P-p0))^1/2;
    dmdt = mdot;
else
    dmdt = 0;
end

