clc
clear

% Exercise 6

S = [[0;1],[1;0],[1;1],[0;1]] %initial vertices

% Case (a)
S1 = [[1;1],[2;1],[2;2],[1;2]] % tranfromed  

syms a11 a12 a21 a22 v1 v2

% transformation matrix defined as symbolic unknown variables
A =  [a11,a12;
     a21,a22]
 v = [v1;v2]
 % for case a it is obvious that v = [1;1] as the [0;0] is
 % transformed into [1,1]. No linear 2x2 map can do that exept 
 % by adding constant biase

 
 %product = A*S +v
 S(:,1) =[];
 S1(:,1) = [];
 product = A*S +[1;1]
 
 %derive system of equations:
eq = product == S1
solution = solve(eq)

%lets print out solved A matrix
A  = [solution.a11,solution.a12;
    solution.a21,solution.a22]
% it is identity matrix



% -----------CASE(b)
clc
clear
syms a11 a12 a21 a22 v1 v2 v

A =  [a11,a12;
     a21,a22]

S = [[0;0],[1;0],[1;1],[0;1]] %initial vertices
S1 = [[2;0],[3;0],[4;1],[3;1]]


% in a similar way we see that [0;0] transformed into [2;0]
% which is only possible by applying constant biase [2;0]
v = [v1;v2]
v1 = 2
v2= 0
subs(v)
S(:,1) =[];
S1(:,1) = [];
product = A*S +v
eq = product == S1

solution = solve(eq)

%lets print out solved A matrix
A  = [solution.a11,solution.a12;
    solution.a21,solution.a22]




 