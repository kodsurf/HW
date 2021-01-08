clc
clear
%Rivet calculation


% Extract reaction forces in joints local coordinate systems
% P1 and P4 -x  y -z
%RP1 = [-24.433; 434.81; 1400.9]  %P1 [x;y;z] components [N] 
%RP2 = [10.351;452.24;1747.5] % P4
%RP3 = [24.434;434.83;-1400.9] % P5
%RP4 = [-10.35;452.22;-1747.5] % P8 

RP1 = [-31.374; 379.15; 1345.9]  %P1 [x;y;z] components [N] 
RP2 = [14.538;455.62;1641.2] % P4
RP3 = [31.32;379.08;-1345.9] % P5
RP4 = [-14.552;455.68;-1641.2] % P8 

% 1) Calculate transferred axial and shear forces for every joint


P1_axial = RP1(1)
P1_sheer = sqrt(RP1(2)*RP1(2)+RP1(3)*RP1(3))

P2_axial = RP2(1)
P2_sheer = sqrt(RP2(2)*RP2(2)+RP2(3)*RP2(3))

P3_axial = RP3(1)
P3_sheer = sqrt(RP3(2)*RP3(2)+RP3(3)*RP3(3))

P4_axial = RP4(1)
P4_sheer = sqrt(RP4(2)*RP4(2)+RP4(3)*RP4(3))


% Calculate reserve factor

% Set parameters

Fbr = 300000000
Fs = 120000000
Ft = 180000000 % 180 MPa
F13 = 30000000 % 30 MPa
Fbr_fastener = 350000000 % 350 MPa
Fa_fastener = 280000000% 280 MPa
Fs_fastener = 600000000 %600 MPa


%Geometrical 
e = 0.0125 %mm
s = e
tl =0.0001 % Laminate thicknes
t_wing = 0.0022 % mm wing thichness
t_centerbox = 0.0018 % mm
t = t_wing + t_centerbox


% Calculations

% Strap strength values
%syms d
d =0.004
%Pbr = Fbr*d*tl % Bearing
%Ps = (2*e-d)*Fs*tl % Shear out
%Pt = (2*s-d)*Ft*tl % Net tention
%Ppt = F13*d*pi*tl % Pull-Through

Pbr = Fbr*d*t_centerbox % Bearing
Ps = (2*e-d)*Fs*t_centerbox % Shear out
Pt = (2*s-d)*Ft*t_centerbox % Net tention
Ppt = F13*d*pi*t_centerbox % Pull-Through


%Reserve factor for Fbr
RFbr1 = Pbr/P1_sheer 
RFbr2 = Pbr/P2_sheer
RFbr3 = Pbr/P3_sheer
RFbr4 = Pbr/P4_sheer

% Reserve factors for RFs

RFs1 = Ps/P1_axial
RFs2 = Ps/P2_axial
RFs3 = Ps/P3_axial
RFs4 = Ps/P4_axial

% Reserve factors for RFt

RFt1 = Pt/P1_axial
RFt2 = Pt/P2_axial
RFt3 = Pt/P3_axial
RFt4 = Pt/P4_axial

%Reserve factors for RFpt 

RFpt1 = Ppt/P1_axial
RFpt2 = Ppt/P2_axial
RFpt3 = Ppt/P3_axial
RFpt4 = Ppt/P4_axial


% Fastener Strenght values and Reserve factors

Pbr_fastener = Fbr_fastener*d*t % Bearing
Pa_fastener = Fa_fastener*pi*d*d/4 % Tensile strength
Ps_fastener = Fs_fastener*pi*d*d/4


% Reserve factors for Fbr_fastener

RFbr_fastener1 = Pbr_fastener/P1_sheer
RFbr_fastener2 = Pbr_fastener/P2_sheer
RFbr_fastener3 = Pbr_fastener/P3_sheer
RFbr_fastener4 = Pbr_fastener/P4_sheer


% Reserve factors for Fa_fastener
RFa_fastener1 = Pa_fastener/ P1_axial
RFa_fastener2 = Pa_fastener/ P2_axial
RFa_fastener3 = Pa_fastener/ P3_axial
RFa_fastener4 = Pa_fastener/ P4_axial

% Reserve factors for Fs_fastener

RFs_fastener1 = Ps_fastener/P1_sheer
RFs_fastener2 = Ps_fastener/P2_sheer
RFs_fastener3 = Ps_fastener/P3_sheer
RFs_fastener4 = Ps_fastener/P4_sheer

% Shear and tension interaction

RFinteraction1 = 1/sqrt( (1/RFs_fastener1)* (1/RFs_fastener1)+ (1/RFa_fastener1) *(1/RFa_fastener1))
RFinteraction2 = 1/sqrt( (1/RFs_fastener2)* (1/RFs_fastener2)+ (1/RFa_fastener2) *(1/RFa_fastener2))
RFinteraction3 = 1/sqrt( (1/RFs_fastener3)* (1/RFs_fastener3)+ (1/RFa_fastener3) *(1/RFa_fastener3))
RFinteraction4 = 1/sqrt( (1/RFs_fastener4)* (1/RFs_fastener4)+ (1/RFa_fastener4) *(1/RFa_fastener4))


%solve system of linear equations


Reserve_Factors = [RFbr1;RFbr2;RFbr3;RFbr4;
    RFs1;RFs2;RFs3;RFs4;
    RFt1;RFt2;RFt3;RFt4;
    RFpt1;RFpt2;RFpt3;RFpt4;
    RFbr_fastener1;RFbr_fastener2;RFbr_fastener3;RFbr_fastener4;
    RFa_fastener1;RFa_fastener2;RFa_fastener3;RFa_fastener4;
    RFs_fastener1;RFs_fastener2;RFs_fastener3;RFs_fastener4;
    RFinteraction1;RFinteraction2;RFinteraction3;RFinteraction4]


RFbr = [RFbr1;RFbr2;RFbr3;RFbr4]
RFS = [RFs1;RFs2;RFs3;RFs4]
RFt = [RFt1;RFt2;RFt3;RFt4]
RFpt = [RFpt1;RFpt2;RFpt3;RFpt4]
RFbr_fastener = [RFbr_fastener1;RFbr_fastener2;RFbr_fastener3;RFbr_fastener4]
RFa_fastener = [RFa_fastener1;RFa_fastener2;RFa_fastener3;RFa_fastener4]
RFs_fastener = [RFs_fastener1;RFs_fastener2;RFs_fastener3;RFs_fastener4]
RFinteraction = [RFinteraction1;RFinteraction2;RFinteraction3;RFinteraction4]
logical = abs(Reserve_Factors) > 1
%equation = Reserve_Factors >= 0
%solution = solve(equation,d)


