function [ISP, Tpeak, t] = ISPcalc(force)
    %ISP calc calculates:
    %**the Isp of the bottle rocket with the isp equation outlined in the lab document
    %**the peak thrust of the rocket for the provided test
    %**the total thrusting time 
   
    g       = 9.81;   %m/s^2
    mprop   = 1;      %kg
    
    [R,~]   = size(force);
    t       = R/1652;           %seconds
    x       = linspace(0,t,length(force(:,3)));
    
    intF    = trapz(x,force(:,3));
    
    ISP     = intF/(mprop*g);   %seconds
    Tpeak   = max(force(:,3));  %N
    
end

