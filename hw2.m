clc
clear
% 2.1
syms streamline(r,theta)
syms potential(r,theta)

streamline(r,theta) =  2*r*r*sin(2*theta)



Vr = diff(streamline,theta)*1/r
Vtheta = -diff(streamline,r)

r = 1
theta = pi/4

% a)
V = [subs(Vr); subs(Vtheta)]



% b)     ------------------------------------
%clc
%clear
syms potential(r)
syms theta
syms r

% Deriving using Vr
% Vr = d phi/ d r   = 4*r*cos(2*theta)
% d phi = 4*r*cos(2*theta) *dr 
potential(r) = int(Vr,r)


% Deriving using Vtheta 
% Vtheta = 1/r * d Phi/d theta = -4*r^2*cos(2*theta)

potential2(theta) = int( Vtheta*r, theta )

% which is the same result for potential function. 

% c) 
% now plot stream lines and potential lines

% generate vector field based on stream line function
theta = linspace(0,2*pi,100);
r = linspace(0,6,100);
[theta,r] = meshgrid(theta,r);
VthetaMesh = -4*r*sin(2*theta);
VrMesh = 4*r*cos(2*theta);

%convert to cartesian
[X,Y] = pol2cart(theta,r);
[U,V] = pol2cart(VthetaMesh,VrMesh);

% plot vector field
figure(1)
hold on
title('Velocity vector field (Cartesian)')
axis equal
quiver(X,Y,U,V)
hold off



% plot pottential field
theta_space = linspace(0,2*pi,360);
r_space = linspace(0,6,360);



%potential_field = r_space.*r_space.*2.*(4.*cos(theta_space).^2 - 2)

[x_space,y_space] = pol2cart(theta_space,r_space);
[X,Y] = meshgrid(x_space,y_space);

% lets convert polar function of potential : r^2*(4*cos(theta)^2 - 2)
% into cartesian form:
%potential_field_cartesian = r^2*(4*cos(theta)^2 - 2)
potential_field_cartesian = (X.*X+Y.*Y).*(4.*cos(atan(Y./X)).^2 - 2);


% Potential field as a surface
figure(2)
title('Potential field as surface')
hold on
surf(X, Y, potential_field_cartesian);
grid on
hold off

figure(3)
hold on
title('Potential field')
% potential field as a counor lines
contour(X,Y,potential_field_cartesian,50)
hold off


% Now lets plot stream lines
theta_space = linspace(0,2*pi,360);
r_space = linspace(0,6,360);
[x_space,y_space] = pol2cart(theta_space,r_space);
[X,Y] = meshgrid(x_space,y_space);

%Again convert streamline function from polar to cartesian
% streamline(r,theta) =  2*r*r*sin(2*theta)
streamline_cartesian = (X.*X+Y.*Y).*sin(2.*atan(Y./X));

figure(4)
hold on
title('Stream lines')
contour(X,Y,streamline_cartesian,50)
hold off



%looks good, lets try combining both plots into 1
figure(5)
hold on
title('combined')
contour(X,Y,streamline_cartesian,50)
contour(X,Y,potential_field_cartesian,50)

hold off
% we can see that ploted countors of stream and potential lines are
% perpendiculat to each other. Therefore we can conclude that results are
% valid


%differential equations tutorial
syms streamlines(theta,r)
syms potential(r)
syms theta
streamline(r,theta) =  2*r*r*sin(2*theta)

Vr = diff(streamline,theta)*1/r
Vtheta = -diff(streamline,r)

equation = diff(potential,r) == 4*r*cos(2*theta)
solution_potential = dsolve(equation)