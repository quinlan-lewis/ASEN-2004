%for design iteration of your glider
%spits out v, L/D, and W/S
clear
clc
close all

%% FILL THIS OUT

B = 1;%span [m] max is 1
B2 = .6;
B = B+B2;
SPlan = .1129;%[m^2]
SPlan2 = .0403;
SPlan = SPlan+SPlan2;
ChordRoot = .25; %[m].1270 minimum, 
ChordTip = .08; %[m] .0762 minimum
ChordRoot2 = .1;
ChordTip2 = .1;
ChordRoot = ChordRoot+ChordRoot2;
ChordTip = ChordTip+ChordTip2;
%Using clark Y-14 airfoil
PThick = .117; %Percent thickness
ClMinD = .1; %Minimum drag Cl form airfoil data
MinCd = .029; %minimum Cd from 2d airfoil

%% Preliminary calculations and constants
%Constants
rho = 1.0581; %[kg/m^3]
CeffSkinFric = .0019;
g = 9.8;

%calcs
AR = B^2/SPlan;
AR2 = B2^2/SPlan2;
SWet = WettedAreaQL(PThick, ChordRoot, ChordTip, SPlan);
% SWet2 = WettedAreaQL(PThick, ChordRoot2, ChordTip2, SPlan2);
Weight = WeightCalcQL(B, ChordTip, ChordRoot, PThick); %[N]
% Weight2 = WeightCalcQL(B2, ChordTip2, ChordRoot2, PThick);

%aerodynamic parameters
e0 = 1/((1/(0.99*SWet))+(0.38*pi*MinCd*AR));
% e02 = 1/((1/(0.99*SWet2))+(0.38*pi*MinCd*AR2));
k = 1/(pi*e0*AR);
% k2 = 1/(pi*e02*AR2);

%% Calcuations

CDmin = CeffSkinFric*(SWet/SPlan); 
CD0 = CDmin + k*ClMinD;
CL = sqrt(CD0/k);
CD = CD0 +k*CL^2;

% CDmin2 = CeffSkinFric*(SWet2/SPlan2); 
% CD02 = CDmin2 + k2*ClMinD;
% CL2 = sqrt(CD02/k2);
% CD2 = CD02 +k2*CL2^2;

%returns
v = sqrt(2*Weight/(rho*SPlan*CL))
LoverD = CL/CD
WingLoading = Weight/SPlan

%% Making the Pretty Plot
V=(1:.1:20);

WtoSref=CL.*(.5*rho.*V.^2);
plot(V,WtoSref,'black','linewidth',2);
hold on
plot(v,WingLoading,'r*','MarkerSize',10,'linewidth',4)

vrangex=[10 15 15 10];
vrangey=[0 0 100 100];
patch(vrangex,vrangey,'green','FaceAlpha',.3,'EdgeColor','none')

xlabel('Velocity')
ylabel('W/S')



