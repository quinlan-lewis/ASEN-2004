function Swet_T = WettedAreaQL(tcr, cRoot, cTip, Splf)
%ASEN 2004 Lab 1
% Code for estimating the surface area of the Tempest UA

%main wing wetted

%All DIMENTIONS ARE IN METERS!
%assumes no twist

%Do not touch
tau = tcr/tcr;

%Do not touch
gamma = cTip/cRoot;


%Calculations
Swet_wing = 2*Splf*(1+0.25*tcr*(1+tau*gamma)/(1+gamma));

%horizontal tail wetted
%assume airfoil shape for horizontal plane
gammaH = 1; %ignoring curve
SplfH = .25*.03; %slightly lowered to accoutn for curve
Swet_horiz = 2*SplfH*(1+0.25*tcr*(1+tau*gammaH)/(1+gammaH));

%vertical tail 
%assume flat plate
%assume rectangular
Swet_vert = .02*.04;

%split the body into sections
%nose cone - .15 m long
%center - chord +.05
%tail cone - 1.56-.43-.25 m long = .88 m 

% r = .09/2;
% 
% %center 
% hC = cRoot+.05;%center height
% Swet_cyl = 2*pi*r*hC;
% 
% %nose
% hN = .15;
% Swet_nose = pi*r*sqrt(hN^2+r^2);
% 
% %tail
% hT = .4;
% Swet_tail = pi*r*sqrt(hT^2+r^2);

Swet_fuleslage = (.03*sqrt(.02^2+.01^2)) + (.4*.03);

Swet_T = Swet_fuleslage + Swet_vert + Swet_horiz + Swet_wing;

end
