clc
clear


% Plot power required and power avaliable 


% initial parameters
m0 = 435 % kg
S = 3.51 % m"2 %wing planform area
g = 9.8 % m/s"2
rho0 = 1.225 %kg/m"3    SEA LEVEL STANDART CONDITIONS
p0 = 101.325 * 1000 % Pa
W= m0*g   % At least 4.263 kN


rho = rho0
h=0
[Temp, a, Pressure, rho] = atmosisa(h) % Model for standart atmospheric conditions
h = densityalt(rho)%  converts from rho to altitude in m




% Power avaliable  Pa = Tr*V  ,where
% Pa - power avaliable
%Tr - thrust required
% V - velocity of free stream air

% From previous homework we assumed Thrust avaliable at sea level is
% constant with respect to velocity


V =20:1:170
maxThrust = 900 % N
Tr_avaliable(1:size(V,2)) = maxThrust;

Pa = Tr_avaliable.*V

figure(1)
hold on
title('Power Avaliable [W] ')
grid on
grid minor
xlabel('V [m/s]') 
ylabel('Pa [W]')
plot(V,Pa)
hold off


% Power required  Pr = D*V,
% Where Pr -power required
% D(v) drag force depending on velocity


%Load sea level drag force from previous homework
load("Drag_v_sea_level.mat") %Loaded as D_V
Pr = D_V.*V

figure(2)
hold on
title('Power required [W] ')
grid on
grid minor
xlabel('V [m/s]') 
ylabel('Pr [W]')
plot(V,Pr)
hold off


%Plot combined

figure(3)
hold on
title('Power required and avaliable [W] ')
grid on
grid minor
xlabel('V [m/s]') 
ylabel('P[W]')
plot(V,Pr)
plot(V,Pa)
hold off




% Calculate exess power

Power_exess = Pa-Pr
figure(4)
hold on
title('Exess power [W] ')
grid on
grid minor
xlabel('V [m/s]') 
ylabel('P[W]')
plot(V,Power_exess)
plot(V,Pr)
plot(V,Pa)
hold off



% Calculate possible rate of climb depending on velocity 

RC = Power_exess./W

%Plot rate of climb with respect to velocity
figure(5)
hold on
title('Rate of Climb [m/s] ')
grid on
grid minor
xlabel('V [m/s]') 
ylabel('Rate of Climb [m/s]')
plot(V,RC)

hold off



% Determine maximum climb rate graphically 

[MaxClimbRate,MaxClimbRateIndex] = max(RC)

% Check velocity at which maximum Climb rate is possible

V_max_climb_rate = V(MaxClimbRateIndex) 
% 87 m/s


% Now calculate velocity at which maximum climb rate is possible
% analytically and compare results


% from previous homework we know 
S = 3.51 % m"2 %wing planform area
rho0 = 1.225 %kg/m"3    SEA LEVEL STANDART CONDITIONS
Cd_0 = 0.02 % zero lift drag coefficient
rho0 = 1.225 %kg/m"3    SEA LEVEL STANDART CONDITIONS
K = 0.062 % Drag slope
T = maxThrust % 900 N  
% Thrust avaliable is constant depending on velocity

syms ClimbRate(Vel)

ClimbRate(Vel) =Vel*( T/W- 1/2*rho0*Vel*Vel *(S/W) * Cd_0-W*2*K*1/(S*rho0*Vel*Vel) )

dClimbRate_dVel = diff(ClimbRate,Vel)
equation = dClimbRate_dVel == 0
solution = solve(equation)

%Velocity_of_maximum_climb_rate = solution % -86.7042

check = subs(equation,86.7042)  % result is correct 

% Now lets substitute found value of velocity at which maximum climb rate
% occures into equation, to find value of maximum climb rate analytically

Max_Climb_Rate_analytical = ClimbRate(86.7042) %10.3127

% Answer - analytically determined value of velocity at which maximum climb
% rate occures is 86.7042 m/s


% For extra verification lets plot analytical version of Climb Rate
% equation and compare its plot with what we previously obtained

RC_analytical = V.*( T/W- 1/2*rho0.*V.*V .*(S/W) .* Cd_0-W.*2.*K.*1./(S.*rho0.*V.*V) )

% Note that we were calculation possible rate of climb 
%GIVEN THAT CLIMB ANGLE = 0 thus , cos(climb_angle) = 1




%Plot analytical rate of climb
figure(6)
hold on
title('Rate of Climb [m/s] ')
grid on
grid minor
xlabel('V [m/s]') 
ylabel('Rate of Climb [m/s]')
plot(V,RC)
plot(V,RC_analytical)
%plot point of maximum climb rate 
scatter(V(MaxClimbRateIndex),RC(MaxClimbRateIndex),'filled')
text(V(MaxClimbRateIndex),RC(MaxClimbRateIndex)+1,"Maximum climb rate : 10.31 m/s")
text(V(MaxClimbRateIndex),RC(MaxClimbRateIndex)+ 3,"V : 86.70 m/s")

hold off
% Conclusion - results of plots are identical 





% Calculate maximum climb angle


% Lets derive horizontal velocity from knowing total velocity V and
% vertical velocity RC

% RC = V*sin(thetea)
% where theta - climb angle.
% First lets derive theta for every value in RC 

theta_rad  = asin(RC./V)
% convert to degres for more intuitive visualisation
theta_deg = rad2deg(theta_rad)

% lets plot theta with respect to climb rate for debug
figure(7)
hold on
title('DEBUG: Climb angle with respect to vertical velocity ')
grid on
grid minor
xlabel('Theta [deg]') 
ylabel('Rate of Climb [m/s]')
plot(theta_deg,RC)
hold off



% Lets derive horizontal compomnent of velocity Vx
% Note that RC is Vy
Vy = RC
Vx = V.* cos(theta_rad)


% Plot vertical velocity with respect to horizontal velocity

figure(8)
%axis equal
hold on
title('Vy with respect to Vx ')
grid on
grid minor
xlabel('Vx [m/s]') 
ylabel('Vy (Rate of Climb) [m/s]')
plot(Vx,Vy)
hold off


% Based on this plot we can derive maximum climb angle and plot it to the
% same graph

% We already know theta_deg array and can find maximum angle from there
[maxTheta,maxThetaindex] = max(theta_deg)
%velocity at maximum climb angle
V_max_climb_angle = V(maxThetaindex)  % 59 m/s
% lets verify that graphically

% climb angle is  atan(Vy,Vx) 
% and the slope is Vy/Vx
% Thus climb angle can be visualized simply as a line starning from origin
% and passing through point [Vx,Vy]


% Calculate climb angle at which maximum climb rate occures
max_climb_rate_theta_rad = atan(Vy(MaxClimbRateIndex)/Vx(MaxClimbRateIndex))
max_climb_rate_theta_deg = rad2deg(max_climb_rate_theta_rad)

% Plot angles at which maximum climb rate occures
% Plot maximum climb angle 
maxClimbRateX =[0,Vx(MaxClimbRateIndex)]
maxClimbRateY =[0,Vy(MaxClimbRateIndex)]

maxClimbAngleX =[0,Vx(maxThetaindex)]
maxClimbAngleY =[0,Vy(maxThetaindex)]

figure(9)
%axis equal
hold on
title('Angle at which maximum climb rate occures, and maximum climb rate angle ')
grid on
grid minor
xlabel('Vx [m/s]') 
ylabel('Vy (Rate of Climb) [m/s]')
plot(Vx,Vy)
plot(maxClimbRateX,maxClimbRateY) % Max Climb Rate angle slope
plot(maxClimbAngleX,maxClimbAngleY) % Max climb angle slope
scatter(Vx(MaxClimbRateIndex),Vy(MaxClimbRateIndex),'filled')
scatter(Vx(maxThetaindex),Vy(maxThetaindex),'filled')
text(Vx(MaxClimbRateIndex)-30,Vy(MaxClimbRateIndex)+2,"Angle at which of maximum climb rate occures 6.80 DEG ")
text(Vx(maxThetaindex)+10,Vy(maxThetaindex),"<-- Maximum possible climb angle 8.08 DEG ")
hold off




%%%%--------- ÑOMPUTE MAXIMUM CLIMB ANGLE ANALYTICALLY


% Use equation RC = V*sin(theta) = exess_power/W

syms theta_func(Vel)
% we know that maximum climb angle is expected to be  8 deg. 
%Thus cos(8deg) = 0.99 , therefore let assume that cos(theta)*cos(theta)
%term in equation is = 1 as we are dealing with small angles only
theta_func(Vel) = asin (  T/W- 1/2*rho0*Vel*Vel *(S/W) * Cd_0-W*2*K*1/(S*rho0*Vel*Vel) )


derivative = diff(theta_func,Vel)
equation = derivative == 0
solution = solve(equation)

max_climb_angle_Velocity_analutically = solution(2) % 59.0869 m/s

% Now lets substitute this value to the function and see what is the
% maximum theta angle calculated analytically
max_theta_analytical_rad = theta_func(max_climb_angle_Velocity_analutically) %0.1412 rad
max_theta_analytical_deg = rad2deg(0.1412) % 8.0902 deg

% exactly the value that we have expected.


%%%% Plot maximum takeoff weight depending on velocity
% We have a graph of power exess depending on velocity V
%Power_exess

% We know that rate of climb RC =  exess_power/W
%RC = Power_exess./W