clc
clear
% Physics of flight hw 4


% Setup drag polar equation

%Cd = 0.02 + 0.062Cl*Cl

Cl = linspace(-0.5,1.4,100)
Cd = 0.02+0.062*Cl.*Cl
figure(1)
title('Drag polar')
%axis equal
plot(Cd,Cl)
title('Drag polar')
grid on
grid minor
xlabel('Cd') 
ylabel('Cl')
%axis equal





% Required thrust at sea level

%take of mass m0 = 435 kg

% In order to take off at sea level lift force L should be greater then
% weight

% L > m*g


% initial parameters
m0 = 435 % kg
S = 3.51 % m"2 %wing planform area
g = 9.8 % m/s"2
rho0 = 1.225 %kg/m"3    SEA LEVEL STANDART CONDITIONS
p0 = 101.325 * 1000 % Pa


L = m0*g   % At least 4.263 kN

% Lets plot graph of Cl(V) required to produce L=4.263 kN at velocity V
% Cl(V) = L /(q*S)
% q = 1/2*rho*V*V

V =20:1:170

q  = 1/2*(rho0*V.*V)
Cl_V = L ./(q.*S)



% Now lets make another plot Cl required to produce constant lift depending
% on velocity V
% This graph would represent minimum takeoff velocity depending on Cl 
%Cl(V) 
%plot(V,Cl_V)
figure(2)
%axis equal
plot(V,Cl_V)
title('Cl required to produce 4.263 kN lift depending on velocity V')
grid on
grid minor
xlabel('V m/s') 
ylabel('Cl')
%axis equal

% Now lets compute drag coefficient Cd(V) which would be produced when flying
% with velocity V and lift coefficient Cl by known relation

%Cd = 0.02 + 0.062Cl*Cl

Cd_V = 0.02+0.062*Cl_V.*Cl_V

% Lets again plot Cd(V) 

figure(3)
%axis equal
plot(V,Cd_V)
title('Cd from producing lift L= 4.263 kN depending on velocity V')
grid on
grid minor
xlabel('V m/s') 
ylabel('Cd')



% Now lets recalculate drag force D depending on velocity V  

D_V = 1/2*rho0*S*Cd_V.*V.*V

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D0_V = 0.02*q*S
D_lift_V = Cl_V.*q*S*0.062
totalD = D0_V+D_lift_V
totalD = 0.02*q*S+2*0.062*S./(rho0*V.*V) *(L*L/S*S)
% and also plot it here
figure(4)
hold on
plot(V,totalD)
title('Total Drag ')
xlabel('V m/s')
ylabel('D [N]')
grid on
grid minor
hold off

% Plot

figure(5)
%axis equal
plot(V,D_V)
title('D(V) Sea level ')
grid on
grid minor
xlabel('V m/s') 
ylabel('D [N]')


maxThrust = 900 % N
maxThrustArray = V;
maxThrustArray(:) = maxThrust
hold on
plot(V,maxThrustArray)

% Required thrust to sustain lift force = 4.263 kN is equal to drag D(V)


% Knowing that minumum required thrust (Drag) occures at maximum L/D ratio
% Lets find minimum drag, and velocity at which it occures

[minD,minIndex] = min(D_V) 
minimumVelocity = V(minIndex)
% Minimun drag = 300N at velocity 59 m/s

% Knowing that maximum avaliable thrust is 900N 
% Lets find maximum velocity that we can achieve at sea level

[minIndex] = find(D_V<905 & D_V>895)
maximumVelocity = V(minIndex(1))
txtString = sprintf("<-max V=%0.2fm/s",round(maximumVelocity,2),round(maxThrust,2));
text((maximumVelocity+5),maxThrust,txtString);
scatter(maximumVelocity,maxThrust,'filled')

scatter(minimumVelocity,minD,'filled')

%txtString = sprintf("<--minimum drag:%0.5fN at V=%0.5f m/s",minD,minimumVelocity);
txtString = sprintf("<--Min drag:%0.1fN at V=%0.1f m/s",round(minD,2),round(minimumVelocity,2));
text((minimumVelocity+10),minD,txtString);

txtString = sprintf("Avaliable thrust: %0.2fN",maxThrust);
text((size(V,2)/2-30),(maxThrust+55),txtString);

hold off



%timestring = sprintf("time passed = %0.5f",total_t);
%text(x_end/2,Vx_x0/2,timestring)

% This was graphical way
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lets try analytically
% For that refere to equations :
%D_V = 1/2*rho0*S*Cd_V.*V.*V
% Cd_V = 0.02+0.062*Cl_V.*Cl_V

%q  = 1/2*(rho0*V.*V)
%Cl_V = L ./(q.*S)

% Derive D(V) function 





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% Lets calculate avaliable thrust depending on altitude

h = 0:100:4000


% Ta/ Ta0 = rho/rho0
% D(rho) = D0*rho/rho0
D0 = maxThrust %N sea level thrust

[Temp, a, Pressure, rho] = atmosisa(h) % Model for standart atmospheric conditions
D_rho= D0*rho./rho0
h = densityalt(rho)%  converts from rho to altitude in m

% Plot thrust required with respect to altitude h in m

figure(6)


plot(h,D_rho)
grid on
grid minor
xlabel('altitude h [m]') 
ylabel('D [N]')

% Calculate avaliable thrust at altitude of 10 000ft
h1 = 3048 % m 
[i, ii, iii, rho1] = atmosisa(h1)
Thrust_h1= D0*rho1/rho0 % 664 N

% Lets plot this  constant line at the same plot 
D_constant(1:size(h,2)) = Thrust_h1
hold on
plot(h,D_constant)


% Now lets recalculate D(V) at altitude of 10 000ft
h1 = 3048
g = 9.8
L = m0*g
V =20:1:170 %  m/s
q  = 1/2*(rho1*V.*V)
Cl_V = L ./(q.*S)
Cd_V = 0.02+0.062*Cl_V.*Cl_V
D_V = 1/2*rho1*S*Cd_V.*V.*V

% Plot drag force ( thrust required) at 3000m  10 000ft

% Plot

figure(7)
%axis equal
plot(V,D_V)
title('D(V)  altitude = 3000m')
grid on
grid minor
xlabel('V m/s') 
ylabel('D [N]')


maxThrust = Thrust_h1
maxThrustArray = V;
maxThrustArray(:) = maxThrust
hold on
plot(V,maxThrustArray)


[minIndex] = find(D_V<Thrust_h1+5 & D_V>Thrust_h1-5)
maximumVelocity = V(minIndex(1))
scatter(maximumVelocity,Thrust_h1,'filled')
txtString = sprintf("<-max V=%0.2fm/s",round(maximumVelocity,2),round(Thrust_h1,2));
text((maximumVelocity+10),Thrust_h1,txtString);

[minD,minIndex] = min(D_V) 
minimumVelocity = V(minIndex)
scatter(minimumVelocity,minD,'filled')

txtString = sprintf("<--Min drag:%0.1fN at V=%0.1f m/s",round(minD,2),round(minimumVelocity,2));
text((minimumVelocity+10),minD,txtString);

txtString = sprintf("Avaliable thrust: %0.2fN",maxThrust);
text((size(V,2)/2-30),(maxThrust+55),txtString);

hold off