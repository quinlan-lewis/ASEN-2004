function [weightTotal] = WeightCalcQL(B, ChordTip, ChordRoot, PThick)

%Weight of total glider estimate

%density of foam
rho_foam = 26; %km/m^3
rho_balsa = 160; %kg/m^3
weight_camera = 0.16;%kg

%Volumes for different parts of aircraft.
%V_fuselage = pi*0.045^2*0.15/3 + (pi*0.045^2*(ChordRoot + 0.05))+ (pi*0.045^2*0.40/3);
V_fuselage = .03*((.02*.01*.5)+(.1*.04)+(.06*.8)); %meters cubed

V_wing = (ChordTip + ChordRoot)/2*B*(PThick*(ChordTip + ChordRoot)/4);

V_tail = (.03 + .02)/2*.25*(PThick*(.03 + .02)/4) + (.02*.04*.03);

weightTotal = ((V_fuselage+V_wing+V_tail)*rho_foam + weight_camera)*9.81;
Price = ((V_fuselage+V_wing+V_tail)*rho_foam) * 1000

end