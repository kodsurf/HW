clc
clear
% exe3

syms p(x) a1 a2 a3 x 

%space of linear polynomials of a degree 2 is defined by
p(x) = a1+a2*x+a3*x*x

% lets assume that we have a basis :
pb1 = [1;0;0]
pb2 = [0;x;0]
pb3 = [0;0;x*x]

% thus polenomial can be represented as a vector p
p1 = [a1;a2;a3]


% (b) : 

% the transformation is defined as a derivative of polynomial
% lets try to figure out vector valued function that transforms a1 a2 a3

syms B(a1,a2,a3,x)

%B(a1,a2,a3,x) = [0;a2;2*x*a3]
B(a1,a2,a3) = [a2;2*a3;0]

%lets compute jacobian 
B_j = jacobian(B)
%Which would be our transformation matrix

% lets test it by creating abrivitory p vector
p_test = [1;2;3] % it looks like  p = 1+2x+3x*x
% derivative shoul be p'= 2+6x
%lets see what we get if we multiply it by a obtained jacobian matrix

result = B_j*p_test % and it is [2;6;0] as expected
%That is indeed work

%lets check if this transformation is linear in the same way as in exe 2
is_linear = check_if_transform_is_linear(B) % if 1 function is linear

% (a)
%similarly lets compose function for tranformation
syms A(a1,a2,a3)
A(a1,a2,a3) = [a1;a2;a1+2*a2-a3]
%lets see tha transformation matrix
A_jac = jacobian(A)
%Lets check if transformation is linear
is_linear = check_if_transform_is_linear(A)
% function returns 1 - transformation is linear



% (c)
syms C(a1,a2,a3)

% figure out transformation matrix
C(a1,a2,a3) = [0;a1*a2;2*a2+a1]

% calculate jacobian matrix
C_jac = jacobian(C)
% Jacobian matrix has non constant values - therefore transformation is not
% linear

%lets confirm that by checking analytically 
is_linear = check_if_transform_is_linear(C)
%function returns 0 
% transformation is indeed non linear


function ret = check_if_transform_is_linear(func)
% Now try analytically prove L(a*u+v) = a*L(u)+L(v)  (knowing that it is
% linear)
syms ux uy uz
syms vx vy vz
u = [ux;uy;uz]
v = [vx;vy;vz]
a=5 % set a to some value
% a*u+v =:
% 5*ux + vx
 %5*uy + vy
 %5*uz + vz
left_side_of_equation = func( 5*ux + vx,5*uy + vy,5*uz + vz)
right_side_of_equation = func(ux,uy,uz)*5 + func(vx,vy,vz)
result = left_side_of_equation - right_side_of_equation

if result == [0;0;0]
    ret=1 %linear
else
    ret=0 %non linear
end


end
