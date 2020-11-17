clc
clear
%exercise 2

syms x y z A(x,y) B(x,y,z) C(x,y,z)

%(a)
%In this case we are given with vector valued functions
A(x,y) = [2*x;-x+3*y;0]

% if Jacobian matrix is constant - therefore transformation might be
% linear
A_map =jacobian(A)

% Exersize states that it it is a transformation from R2 to R3
r_a = rank(A_map)
% However rank of jacobian matrix = 2, therefore it is a transformation
% onto 2D plane in 3D space. 
% which is logical because z component of vector valued function A(x,y) is
% constant 0, thus points from R2 would be mapped onto R3 plane with z=0


%(b)
B(x,y,z) = [x+2*z;
           y-z+2;
            3*z]

        
%All the values of jacobian matrix are again constants, therefore map is
%might be linear (HOWEVER IT IS NOT)
B_map =jacobian(B)

% We can see that rank of jacobian matrix = 3, therefore it is indeed
% transformation from R3 to R3
r_b = rank(B_map)



% Lastly lets prove analytically that transformation is non
% linear by defination given in lecture notes 
% which is :
% Map is linear if 1) L(a*u+v) = a*L(u)+L(v)
% consider case (b) with equation :
clear x y z
syms x y z
B(x,y,z) = [x+2*z;y-z+2;3*z]

% lets define u and v
syms ux uy uz
syms vx vy vz

u = [ux;uy;uz]
v = [vx;vy;vz]

a=5 % set a to some value 

% Now try analytically prove L(a*u+v) = a*L(u)+L(v)  (knowing that it is
% linear)
% a*u+v =:
% 5*ux + vx
 %5*uy + vy
 %5*uz + vz
left_side_of_equation = B( 5*ux + vx,5*uy + vy,5*uz + vz)



right_side_of_equation = B(ux,uy,uz)*5 + B(vx,vy,vz)


%Now lets sum up left and right side of equtions to see if they are equal
% If they are - this is linear map
% If not - mapping is non linear 

result = left_side_of_equation - right_side_of_equation
% We see that result is non zero vector because equations are not equal. 
% the reason for that is because Y(x,y,z) in B(x,y,z) has constant biase 2
% when we calculate jacobian by taking derivatives. Derivative of a
% constant were 0 and bias of 2 were eliminated from equation. 



% (C)
C(x,y,z) = -x+4*y-3*z/2

% as it is a R3 to R1 transformation that linear map is a vector
C_map = jacobian(C)

% lets apply this transformation for a demo purposes:

%example vector in R3
r3 = [1;2;3]

%transformed using matrix
r3_transformed = C_map*r3

% Now lets transform the same vector analytically: 
% 
x=1
y=2
z=3

r3_transformed_analyticaly =  C(x,y,z)
%from a console command we se that results are equal
is_equal = r3_transformed == r3_transformed_analyticaly


%We can try checking that this is a linear map and jacobian matrix is
%indeed transformation matrix. However I can tell already in advance that
%it is because there are no biases in C(x,y,z)

% lets define u and v
clear ux uy uz vx vy vz
clear x y z left_side_of_equation right_side_of_equation
syms x y z
syms ux uy uz
syms vx vy vz
u = [ux;uy;uz]
v = [vx;vy;vz]
a=5 % set a to some value 

% Now try analytically prove L(a*u+v) = a*L(u)+L(v)  (knowing that it is
% linear)
% a*u+v =:
% 5*ux + vx
 %5*uy + vy
 %5*uz + vz
left_side_of_equation = C( 5*ux + vx,5*uy + vy,5*uz + vz)
right_side_of_equation = C(ux,uy,uz)*5 + C(vx,vy,vz)

result = left_side_of_equation - right_side_of_equation
% Now it is indeed 0 there fore transformation is a linear map