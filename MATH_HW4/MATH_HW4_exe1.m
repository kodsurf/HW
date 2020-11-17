clc
clear
%exersice 1

% GIVEN : set of vectors b and p
b1 = [-2;0]
b2 = [3;2]

p1 = [1;1;0]
p2 = [-2;-2;-2]
p3 = [3;0;2]


% lets compose matrixes out of b and p basic vectors

B = [b1,b2]
P = [p1,p2,p3]

%I prefere to denote vectos are row vectors. But lets transpose to column 
%vectors as in task defination
B = transpose(B)
P = transpose(P)


B_rank = rank(B) % 2 : indeed can be basis for R2
P_rank = rank(P) % 3 : indeed can be basis for R3
% I am not doing gausian ellimination manually to show that this is linear
% independant vectors. As this was in previous HW.  Rank is sufficient to
% demonstrate that vectors can be a basis. 

% lets define vector v analyticaly     v = [1;2] with respect to B
v = 1* b1+ 2*b2   % where b1 b2 are basis vectors and 1 ,2 components with respect to B basis

% now we can already see that v vector with respect to standart basis is
% [4;4]

%and for P it is similar :  Where v and w are vectors with respect to
%standart basis. And both v and w were defined by a components of B and P
%basis
w = 4*p1+1*p2+6*p3 % [20;2;10]

% Now lets assume that v and w are defined with respect to standart BASIS
v_standard = [1;2]
w_standard = [4;1;6]
%and we want to know coordinates of v and w with respect to B and P :v_B,
%v_P

% for that we need to compose FORWARD transformation matrixes from standard
% to B and P basis

% whic is :

B_forward = transpose(B)
P_forward = transpose(P)

%Numerically this forward transformation converts vector FROM B to standart
% FROM P to standart.
% However it is still called forward transformation.
%lets verify :

result = B_forward*[1;2]
% we got [4;4]  with respect to standart basis
% in the same way we obtained [4;4] when calculated analytically 
%v = 1* b1+ 2*b2   

result = P_forward * [4;1;6]
% and again we got [20;2;10] with respect to standart basis as in analytical defination:
%w = 4*p1+1*p2+6*p3 % [20;2;10]



% lets move forward with our task and compute BACKWARD transformation
% matrix

B_backward = inv(B_forward)
P_backward = inv(P_forward)

% whicj is inverse of forward transformation matrix
% numericaly bacward transformation matrix calculates FROM standard to new
% B and P basises

% lets prove this numerically

v_B = B_backward * v_standard %[1;1]
w_P = P_backward* w_standard  %[-3;-2;1]

% this is a coordinates of vectors v and w with respect to B and P basises

% lets try converting v_B= [1;1] and w_P ] [-3;-2;1] BACK into standard
% basis but this time based on analytical vector defination:
v_standart = b1*1 + b2*1 % [1;2]
v_standart = p1*-3 +p2*-2 + p3*1 %[4;1;6]

% and we indeed get the same vector


% Conclution : 
% 1) VECTORS are always defines as linear combination of COMPOTENTS and a
% BASIS vectors
%2) The same vector can be defined with DIFFERENTS compotents depending on
%BASIS which are used to describe the vector. 
