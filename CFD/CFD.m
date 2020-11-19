clc
clear
% CFD 

%Lets consider 1D example of conservation of momentim equation:
% rho - density
% V = [vx;vy;vz] - velocity vector field ( vector valued function)
% in our case V = [vx(x)]

% lets writte down left part of equation:
% rho pd(vx)/pdt + rho * dot(V ,  gradient(vx))  

%devide by rho we obtain
% pd(vx)/pd(t) + dot(vx,gradient(V))

% for 1D case gradient of V is simply partial derivative with respect to x
% pd(vx)/pd(t) + dot(vx, pd(vx)/pd(x))

% and a dot product of two scalar functions is simply multiplication of
% those functions

% pd(vx)/pd(t) + vx* pd(vx)/pd(x)

% equation describes change in momentum which is equal to summ of forces
% ma = f

% however lets assume that there are not external forces acting in this
% example for simplicity

%pd(vx)/pd(t) + vx* pd(vx)/pd(x) = 0

% This equation os called Burgins equation

syms vx(x,t) t x

equation = diff(vx,t) + vx* diff(vx,t) ==0

% Problem statement:

% Straight pipe  l =2m
%  A------------------------x=0----------------------------B--->x

% Set origin of coordinate system somewhere in between A and B

%  Boundary conditions :
% velocities vx(A,t>=0) = 1.2 m/s     vx(B,t>=0) = 0.4 m/s   are given (for
% every moment in time

% Note that there is a discontinuety 
% we set  vx(  A<x <0, t=0) = vx(A) 
% and vx( 0<x < B,t=0) = vx(B)
% velocities of every cell LEFT to origin at TIME t = 0 is equal to
% vx(A)
% velocities of every cell RIGHT to origin at TIME t=0 is equal vx(B)



%Now lets create cells X
real_cells = 100
x0 = -2
x_end = 6
dX = (abs(x0 -x_end))/(real_cells-2)

for i=1:real_cells+2
    X(i) = x0 +(i-3/2)*dX
end


% Create timesteps
tMax = 9
dt = 0.01
time_cells_number = int16( tMax / dt)
t = linspace(0,tMax,time_cells_number);


% Create and assign vx array for every cell in X 
% Note that there is a discontinuety 
% we set  vx(  A<x <0, t=0) = vx(A) 
% and vx( 0<x < B,t=0) = vx(B)
% velocities of every cell LEFT to origin at TIME t = 0 is equal to
% vx(A)
% velocities of every cell RIGHT to origin at TIME t=0 is equal vx(B)
Vx_A = 1.2
Vx_B = 0.4
Vx =[]
logical = X<=0 % if less or equal to 0  cell value =1 else 0
for i= 1:size(X,2)
    if logical(i) ==1 % if we are to the left from origin
        Vx(1,i) = Vx_A
    else % if we are to the right from origin
        Vx(1,i) = Vx_B
        
    end % end if
    
    
end % end for

% Convert analytical derivative into numerical

% pd(vx(x,t) )/ pd(x) = vx(x+1,t) - vx(x,t)/ dx
% where x is index of X cell   ; dx = X(x+1) - X(x)


%pd (vx)/pd(t) = vx(x,t+1) - vx(x,t)/ dt

% Analytical equation : pd (vx)/pd(t) + vx* pd(vx)/ pd(x)=0

% Discreete equation :
% (vx(x,t+1) - vx(x,t) )/ dt   + vx(x,t) * ( (vx(x+1,t) - vx(x,t))/ dx  )=0

% From which we derive vx(x,t+1)
% vx(x,t+1) = vx(x,t) - dt * ( vx(x,t) * (vx(x+1,t) - vx(x,t))/ dx )


% Main loop
figure (1)

for tn =1:time_cells_number-1 % for every timedtep
    dt = t(tn+1) -t(tn);
    %Now set boundaty conditions in ghost cells
    Vx(tn,1) = Vx_A;
    Vx(tn,end) = Vx_B;
    
    for xi = 2:(size(X,2)-1) % for every REAL position cell
        dx = X(xi+1)-X(xi);  
        dx = abs(dx);
        % update Vx at cell with xi index
        Vx(tn+1,xi) = Vx(tn,xi) - dt*Vx(tn,xi) * (  Vx(tn,xi) - Vx(tn,xi-1)  ) /dx ;
        %Vx(tn+1,xi) = Vx(tn,xi) - dt*Vx(tn,xi) * (  Vx(tn,xi+1) - Vx(tn,xi)  ) /dx ;
        
    end % end for 
    
    % plot for every time step (COMMENT to speed up the loop)
    scatter(X,Vx(tn,:));
    text(X(xi),Vx(tn,xi),"we are accelerated")
    text(X(xi-int16(size(X,2)/2)),Vx(tn,xi-int16(size(X,2)/2)),"Oh myyy")
    pause(0.01)
end %end main
clf
text(0.5,0.5,"Damn !!!");
% NOTE :
% update formula for Vx in next time cell is not as we have derived it
% previously

%if we plug originaly derived formula - distribution would not be correct
% Vx(tn+1,xi) = Vx(tn,xi) - dt*Vx(tn,xi) * (  Vx(tn,xi+1) - Vx(tn,xi)  ) /dx ;

% Lets try out with original (wrong) formula just for fun
Vx_wrong(1,:) = Vx(1,:)
tn = 1
xi=1
for tn =1:time_cells_number-1 % for every timedtep
    dt = t(tn+1) -t(tn);
    %Now set boundaty conditions in ghost cells
    Vx_wrong(tn,1) = Vx_A;
    Vx_wrong(tn,end) = Vx_B;
    
    for xi = 2:(size(X,2)-1) % for every REAL position cell
        dx = X(xi+1)-X(xi);  
        dx = abs(dx);
        % update Vx at cell with xi index
        %Vx_wrong(tn+1,xi) = Vx_wrong(tn,xi) - dt*Vx_wrong(tn,xi) * (  Vx_wrong(tn,xi) - Vx_wrong(tn,xi-1)  ) /dx ;
        Vx_wrong(tn+1,xi) = Vx_wrong(tn,xi) - dt*Vx_wrong(tn,xi) * (  Vx_wrong(tn,xi+1) - Vx_wrong(tn,xi)  ) /dx ; % this is origonal formula
        
    end % end for 
    
    
end %end main

% And we indeed see that arround the middle element of Vx_wrong -
% velocities is being accelerated. This is not what we expected.

% This is because we have to calculate numerical derivative by sunstracting
% cell value of Vx to the left instead of to the right

% Change forward numerical derivative into backward and it will work