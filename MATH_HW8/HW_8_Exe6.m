%Math_HW8 Exe6
clc
clear
% exe 6



% Compose vector valued function descriping surface

syms R(r,theta,z) x y z

R(r,theta,z) = [r*cos(theta);r*sin(theta);z]
%compute jacobian
J = jacobian(R)

%tangent vectors are
e_r  = [cos(theta);sin(theta);0]
e_theta = [-r*sin(theta);r*cos(theta);0]
e_z = [0;0;1]


% Compose tangent plane

r = sqrt(2)
theta = pi/4
z = 0
x0 = r*cos(theta)
y0 = r*sin(theta)
z0 = z

normalVector = cross(e_theta,e_z)
vector = [x;y;z] - [x0;y0;z0]


dotProduct = dot(normalVector,vector)
plane_equation = dotProduct ==0
plane_equation = subs(plane_equation)
solution = solve(plane_equation)

%plot tangent plane 


%plot cylinder
r = sqrt(2);
theta = linspace(0,2*pi,30);
x = linspace(-r,r,30);
y = linspace(-r,r,30);
z = linspace (0,2,30);

[THETA,Z] = meshgrid(theta,z);
X = r*cos(THETA);
Y = r*sin(THETA);
[gg,Z] = meshgrid(x,z);
figure(1)
hold on
surf(X,Y,Z)



% for p(r,0,0)
x = linspace(-r,r,30);
y = linspace(-r,r,30);
z = linspace (0,2,30);
[X,Z] = meshgrid(x,z);
Y = X;
Y(:)= r;
surf(X,Y,Z)
hold off



%for p(r,pi/4,1)

r = sqrt(2);
theta = linspace(0,2*pi,30);
x = linspace(-r,r,30);
y = linspace(-r,r,30);
z = linspace (0,2,30);

[THETA,Z] = meshgrid(theta,z);
X = r*cos(THETA);
Y = r*sin(THETA);
[gg,Z] = meshgrid(x,z);
figure(2)
hold on
surf(X,Y,Z)


x = linspace(-r,r,30);
y = linspace(-r,r,30);
z = linspace (0,2,30);
[X,Z] = meshgrid(x,z);
Y = 2-X

surf(X,Y,Z)
