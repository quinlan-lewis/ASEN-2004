function [dvdt] = phase2ODE(t,Cd,At,d_water,P0,gamma,P_amb,Vb,R2,M_btl,P_end,T_end,ls,a_x0,a_z0,vel_x0,vel_z0,Ab,g,Rey,Ve,Mr_0,Mair_i)
    P = P0;
    if P > P_amb
        Mair_0 = Mr_0 - M_btl;
        P = P_end*(Mair_0/Mair_i)^gamma;
        P_star = P*(2/(gamma+1))^(gamma/(gamma-1));
        rho = Mair_0/Vb;
        T = P/(rho*Rey);
        if P_star > P_amb
            Te = (2/(gamma+1))*T;
            Pe = P_star;
            rhoe = Pe/(Rey*Te);
            Ve = sqrt(gamma*R2*Te);
        elseif P_star < P_amb
            Me = sqrt((2*(P/P_amb)^((gamma-1)/gamma)-1)/(gamma-1));
            Te = T/(1+((gamma-1)/2)*Me^2); %for some reason this is a matrix of imaginaries
            %fprintf('Ri = %f\n',Ri)
            %fprintf('P_amb = %f\n',P_amb)
            %fprintf('Te = %f\n',Te)
            rhoe = P_amb/(Rey*Te);
            Pe = P_amb;
            Ve = Me*sqrt(gamma*Rey*Te);
        end
        fprintf('rhoe = %f\n',rhoe)
        fprintf('At = %f\n',At)
        fprintf('Ve = %f\n',Ve)
        mdotair = -Cd*rhoe*At*Ve;
        Thrust = mdotair*Ve + (P_amb - Pe)*At;
        if Thrust < 0
            Thrust = 0;
        end
        %fprintf('%f\n',R(1))
        velx = R2(1);
        velz = R2(2);
        veltot = sqrt(velx^2+velz^2);
        xdir = velx/veltot;
        zdir = velz/veltot;
        xThrust = Thrust*sin(atan(zdir/xdir));
        zThrust = Thrust*cos(atan(zdir/xdir));
        
        Drag = -(rho/2)*(veltot^2)*Cd*Ab;
        xDrag = -(rho/2)*(R2(3)^2)/((cos(atan(zdir/xdir)))^2)*Cd*Ab;
        zDrag = -(rho/2)*(R2(4)^2)/((cos(atan(zdir/xdir)))^2)*Cd*Ab;
        
        mdotrocket = -mdotair;
        
        accx = ((Thrust+Drag)*xdir)/Mr_0;
        accz = ((Thrust+Drag)*zdir)/Mr_0 - g;
        
        velx = R2(1);
        velz = R2(2);
        
        dvdt = [accx,accz,velx,velz,Thrust]';
        
    else
        R2(1) = 0;
        R2(2) = 0;
        R2(3) = 0;
        R2(4) = 0;
        R2(5) = 0;
    end
    
end
