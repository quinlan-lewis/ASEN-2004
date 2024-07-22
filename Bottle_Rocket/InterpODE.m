function [mdot] = InterpODE(~,vq)
    global d_water D_thr At

    d_water = 1000; %density of water kg/m^3
    At = pi*(D_thr/2)^2; %Throat Area
    %Ab = pi*(D_btl/2)^2; %bottle area
    mdot = sqrt(d_water*At*vq);
%     xdotdotdot = vg/mdot;
%     xdotdot = xdotdotdot;
%     xdot = xdotdot;
%     
%     data = [mdot xdotdotdot xdotdot xdot];
%     data = data';
end

