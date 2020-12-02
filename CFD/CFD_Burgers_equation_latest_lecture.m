% UFC Copyright protection ©
clc
clear
% CFD 
%{

%Lets consider 1D example of conservation of momentim equation:
% rho - density
% V = [vx;vy;vz] - velocity vector field ( vector valued function)
% in our case V = [vx(x)]

% lets writte down left part of equation (page 26 CFD):
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
x0 = -pi/2
x_end = pi/2
dX = (abs(x0 -x_end))/(real_cells-2)

for i=1:real_cells+2
    X(i) = x0 +(i-3/2)*dX
end


% Create timesteps
% NOTE :
% In order for solution to be valid time step dt should be selected
% according to Vmax*dt/dx <= 1
% Where Vmax - maximum velocity of fluid particles   in our case it would
% be Vx_A
% dt -time step
% dx - distance between X cells  dx <= |x - X(i)|
tMax = 5
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

% Lets try to generate smooth initial condition settings
currX =0
for i= 1:size(X,2)
    Vx(1,i) = Vx_A*sin(X(i)-pi)
    if Vx(1,i) < Vx_B
        Vx(1,i) = Vx_B
    end
    
    
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
%axis equal
axis ([-30 30 -30 30])
total_t = 0
for tn =1:time_cells_number-1 % for every timedtep
    dt = t(tn+1) -t(tn);
    total_t = total_t+dt
    tt(tn) = total_t;
    %Now set boundaty conditions in ghost cells
    Vx(tn,1) = Vx_A;
    Vx(tn,end) = Vx_B;
    
    for xi = 2:(size(X,2)-1) % for every REAL position cell
        dx = X(xi+1)-X(xi);  
        dx = abs(dx);
        % update Vx at cell with xi index
        %Vx(tn+1,xi) = Vx(tn,xi) - dt*Vx(tn,xi) * (  Vx(tn,xi) - Vx(tn,xi-1)  ) /dx ;
        %Vx(tn+1,xi) = Vx(tn,xi) - dt* (Vx(tn,xi+1) -Vx(tn,xi-1) )/(2*dx)*Vx(tn,xi); % unstable solution with oscilations NO UPWINDING
        %Vx(tn+1,xi) = Vx(tn,xi) - dt/dx * (0.5 * (Vx(tn,xi))^2 - 0.5 * (Vx(tn,xi-1))^2); % correct solution
        

        
        %fVx_p = 0.5 * (Vx(tn,xi))^2;
        %fVx_m = 0.5 * (Vx(tn,xi-1))^2;
        
        if (0.5*(Vx(tn,xi)+Vx(tn,xi+1)) >=0) % flow from left to right fVx_p
            fVx_p = 0.5 * (Vx(tn,xi))^2;
            
            
        else % flow from left to right
            fVx_p = 0.5 * (Vx(tn,xi+1))^2;
            
            
        end
        
        if 0.5*(Vx(tn,xi-1) +Vx(tn,xi)) >= 0 % if flow from left to right fVx_m
            fVx_m = 0.5 * (Vx(tn,xi-1))* 0.5*(Vx(tn,xi-1) +Vx(tn,xi));
            
        else % if flow from right to left
            fVx_m = 0.5 * (Vx(tn,xi+1))* 0.5*(Vx(tn,xi+1) +Vx(tn,xi));
            
        end
        
        
        
        %[fVx_DD,fVx_D] = flux_func(tn,xi,Vx)
        [fVx_p,fVx_m] = flux_func(tn,xi,Vx)
        
        Vx(tn+1,xi) = Vx(tn,xi) - dt/dx * (fVx_p - fVx_m); % correct solution
        
     
        
        
        
        
        
    end % end for 
    
    % plot for every time step (COMMENT to speed up the loop)
    scatter(X,Vx(tn,:));
    text(X(xi),Vx(tn,xi),"we are accelerated")
    text(X(xi-int16(size(X,2)/2)),Vx(tn,xi-int16(size(X,2)/2)),"Oh myyy")
    timestring = sprintf("time passed = %0.5f",total_t);
    text(x_end/2,Vx_A/2,timestring)
    pause(0.01)
end %end main
%clf
%text(0.5,0.5,"Damn !!!");
% NOTE :
% update formula for Vx in next time cell is not as we have derived it
% previously

%if we plug originaly derived formula - distribution would not be correct
% Vx(tn+1,xi) = Vx(tn,xi) - dt*Vx(tn,xi) * (  Vx(tn,xi+1) - Vx(tn,xi)  ) /dx ;





% ∫  ∂  ρ μ  ∆ λ ∇ 

% Lets derive different form of Burgers equation:
%∂Vx/∂t + vx* ∂Vx/∂x = 0  (1.0)  ###
% lets rewrite is as 
% Vx,t + Vx*Vx,x = 0      (1.0)  ###

% Where  Vx,t Vx,t are derivatives with respect to t and x

% Referring to page 35 of CFD Burgers equation can also be defined as:
% Vx,t + f(Vx),x = 0 (1.1) ########3
% Where f(Vx) = 1/2 Vx"2    Flux function 
% We can verify it by taking derivative of ∂f(Vx) /∂x
% ∂f(Vx) /∂x = ∂f(Vx)/ ∂Vx   * ∂Vx/∂x  = 2/2 Vx * ∂Vx/∂x

% Lets take double integral of left and right hand sides of Burgers
% equation

% 0t∫ab∫ (Vx,t + f(Vx),x )dxdt = 0

% Using property - integral of sym is equal to summ of integrals

% 0t∫ab∫ (Vx,t + f(Vx),x )dxdt = 0t∫ab∫Vx,t dxdt +0t∫ab∫f(Vx),x dxdt = 0
% (1.1)###############
% where :
% 0t∫ - definate integral from 0 to t
% ab∫  - definate integral from a to b
% f(Vx) = Vx"2/2 - flux function 
% f(Vx),x - notation for derivative with respect to x  : Note that Vx(x,t)
% x = a , x=b   are boundary conditions (length)

%lets consider second term of right hand side  0t∫ab∫f(Vx),x dxdt 


% we see that f(Vx),x is a derivative of flux function with respect to x,
% therefore integral of the derivative would be the function itself 

%0t∫ab∫f(Vx),x dxdt = 0t∫ ab|f(Vx) dt
%0t∫ab∫f(Vx),x dxdt = 0t∫ [f(Vx(x=b,t)) - f(Vx(x=a,t))] dt

% By substituting f(Vx) = Vx"2 /2 into definate integral we get :
% 0t∫ab∫f(Vx),x dxdt = 0t∫ [(Vx(x=b,t))"2/2 - (Vx(x=a,t))"2/2 ] dt

% Recall boundary conditions : Vx(x=a) = 1.2 m/s  Vx(x=b) = 0.4 m/s 
% By defination velocity is constant at boundary conditions, therefore we
% are integration over constants
% 0t∫ [(Vx(x=b,t))"2/2 - (Vx(x=a,t))"2/2 ] dt = 0t∫ [Vx(b) -Vx(a)]dt =
% = 0t| [Vx(b)"2/2-Vx(a)"2/2]t = [Vx(b)"2/2-Vx(a)"2/2]t

% 0t∫ab∫f(Vx),x dxdt = [Vx(b)"2/2-Vx(a)"2/2]t    (1.2)########

% Note that Vx(b) and Vx(a) are known constants from boundary conditions


% Lets look again on equation :
% 0t∫ab∫ (Vx,t + f(Vx),x )dxdt = 0t∫ab∫ Vx,t dxdt +0t∫ab∫f(Vx),x dxdt
% Second term of right hand side we already derived 
%Focus on right hand side first term 
% 0t∫ab∫ Vx,t dxdt

% We can change the order of integration 
% 0t∫ab∫ Vx,t dxdt = ab∫0t∫ Vx,t dtdx
%As our definate integral limits are constants we can do it withoit
%recalculating intervals 

% https://math.stackexchange.com/questions/1737230/order-of-integration-for-multiple-can-be-easily-swapped-if-limits-are-constants

% ab∫0t∫ Vx,t dtdx = ab∫ [Vx(x,t)-Vx(x,0)]  (1.3) ###########

% Now we again use sum of integrals property
%ab∫ [Vx(x,t)-Vx(x,0)] = ab∫Vx(x,t) - ab∫Vx(x,0)  ## (1.4)


% We are ready to substitute all derived componets into 1.1 equation to get

% ab∫Vx(x,t) - ab∫Vx(x,0) + [Vx(b)"2/2-Vx(a)"2/2]t  = 0  (1.5) ####

% From which we derive 
% ab∫Vx(x,t) = ab∫Vx(x,0) - [Vx(b)"2/2-Vx(a)"2/2]t   (1.6) ####

% However this is continious analytical equation. It wont work applied to
% discrete problem



% We do it by going back to equation 1.1
% 0t∫ab∫Vx,t dxdt +0t∫ab∫f(Vx),x dxdt = 0

% and instead of integration over whole length and time period - we
% integrate over a time cell and cell length

%tn:tn+1∫xi:xi+1∫Vx,t dxdt +tn:tn+1∫tn+1∫f(Vx),x dxdt = 0

% We repeat all the previos steps and end up with discrete equation
% instead a and b we substitude xi-1 and xi+1
% instead 0 and t we substitute tn and tn+1

% finall equation :
% Vx(tn+1,xi) = Vx(tn,xi) - dt/dx * (0.5 * (Vx(tn,xi))^2 - 0.5 * (Vx(tn,xi-1))^2); % correct solution




% Now lets take a look on to right hand side of burgers equation.
% Which are the summ of forces acting on a system
% Previonsly we made a computation assuming that forces are 0
% Lets edit our solution in order to take into account external forces

% ∂Vx/∂t + vx* ∂Vx/∂x = - px/ρ +4/3 * Vx* μ / ρ 

%Right hand side of equations seems to be originated from CFD lecture notes
%page 26:

% ∂(Vx*ρ)/∂t + ∇·(ρVxV) = -∂p/∂x +μ* ∆(Vx)+ (λ+μ) ∂(∇·V)/∂y 
% ∇ · (ρVxV)

%Where :
% V - V(x,y,z) velocity vector
% μ - dynamic viscosity
% λ - "viscosity to relate stresses to the volumetric deformation" CFD
% lecture notes page 24
% p - p(x,y,z) pressure as a function of position

% We are considering liquid incompressible flow therefore :
% ∇ · V = 0 
% Therefore we are neglecting (λ+μ) ∂(∇dotV)/∂y   term of righ hand side
% Incompressible flow also implies that density ρ(x,y,z) is constant,
% therefore we can take it out of equation as a constant value

% ρ*∂(Vx)/∂t + ρ* ∇·(VxV) = -∂p/∂x +μ* ∆(Vx)

% and devide by ρ

% ∂(Vx)/∂t + ∇·(VxV) = (-∂p/∂x)/ρ +μ* ∆(Vx)/ρ

% Lets also transform left hand side of equation to prove that it Burgers
% equations is indeed derived from it

% Focus on a term ∇·(VxV)
% By the chain rule :
% ∇·(VxV) = V* (∇·Vx) + Vx* (∇·V)
% We know that incompressible, therefore (∇·V) = 0

%Moreover we are considering 1D proplem along x axis thus V is function of only x
%possition V(x)  V(x) = Vx

% ∇·Vx divergence if Vx is simply Vx,x   because in 1D Vx(x) is a function
% of only x

% ∂(Vx)/∂t + Vx*(∂Vx/∂x)   = (-∂p/∂x)/ρ +μ* ∆(Vx)/ρ

% Equation defined in 20 lecture is :
%  ∂Vx/∂t + vx* ∂Vx/∂x = - p,x/ρ +4/3 * Vx,xx* μ / ρ 

% Lets integrate right hand side the same way as we have integrated the
% left hand side

% NOTE we are inregrating with a limits of time step and x possition cell
% tn:tn+1∫xi-1:xi+1∫(- p,x/ρ +4/3 * Vx,xx* μ / ρ ) dx dt

% lets start with a term :
% tn:tn+1∫xi-1:xi+1∫- p,x/ρ dx dt

% Straith forward lets derive numerical sheme for thios integral 
% tn:tn+1∫xi-1:xi+1∫- p,x/ρ dx dt = 
% -[ 1/ρ(tn,xi) * p(tn,xi+1)-p(tn,xi-i)/(2 *Δx)]Δx Δt

% Note that here we are using ρ as a function of time and possition, we
% denoted it as a function to use further in a general case of compressible
% flows. However for our numerical solution lets consider that ρ(x,t) is
% constant everythere


%Lets derive numerical scheme for second term 4/3 * Vx,xx* μ / ρ
% tn:tn+1∫xi-1:xi+1∫4/3 * Vx,xx* μ / ρ   dx dt = 
%  4/3*μ/ρ(tn,xi) *[((Vx(tn,xi+1 -2 *Vx(tn,i) + Vx(tn,xi-1))/  Δx"2 ]ΔxΔt

% Finaly combining everything together numerical scheme for updating
% Vx(tn+1,xi)

% Vx(tn+1,xi) = Vx(tn,xi) -dt/dx*(F(xi+1)-F(xi-1) +4/3*μ/ρ(tn,xi)*[((Vx(tn,xi+1 -2 *Vx(tn,i) + Vx(tn,xi-1))/  Δx"2]ΔxΔt
%  -[ 1/ρ(tn,xi) * p(tn,xi+1)-p(tn,xi-1)/(2 *Δx)]Δx Δt



% Now lets clear all the variables and make a fresh code for burgers
% equation considering forces

%}
clc
clear

% ------------ INITIAL PARAMETERS ------------------

% Number of cells
real_cells = 160
ghost_cells = 2 % should be even

% ------ Pipe length 
%x0 = -pi/2
%x_end = pi/2
x0 = -1
x_end =1

% Time of computation and timestep 
tMax = 5
dt = 1e-4


% Boundary conditions for velocities ( Maximum minimum velocity)
Vx_x0 = 1.2
Vx_xend = 0.4


% Viscosity 
const_mu = 5.3

% Density 

const_rho = 1000;
rho_x0 = 1
rho_xend = 0.125

% Specific kinetic energy
e_x0 = 2.5
e_xend = 0

% Gamma
gamma = 1.4


% ---------------------- INITIALIZATION

% Initialize X cells -------------------
dX = (abs(x0 -x_end))/(real_cells);

for i=1:real_cells+ghost_cells
    %X(i) = x0 +(i-3/2)*dX;
    X(i) = x0+(i-(ghost_cells/2+0.5))*dX;
end

% Initialize Time cells 

time_cells_number = int16( tMax / dt)
t = linspace(0,tMax,time_cells_number);


% Initialize Velocity

for i= 1:size(X,2)
    Vx(1,i) = Vx_x0*sin(X(i)-pi)
    if Vx(1,i) < Vx_xend
        Vx(1,i) = Vx_xend
    end
    
    
end % end for

Vx(1,1:int16(size(X,2)/2)) = Vx_x0;
Vx(1,(int16(size(X,2)/2):end)) = Vx_xend;


% Initialize density

for i= 1:size(X,2)
    rho(1,i) = const_rho*sin(X(i)-pi);
end % end for

rho(1,1:int16(size(X,2)/2)) = rho_x0;
rho(1,(int16(size(X,2)/2):end)) = rho_xend;

%Initialize viscosity
muu(1,1:size(X,2)) = const_mu;

% Initialize Specific kinetic energy
e( 1,1:int16(size(X,2)/2) ) = e_x0
e( 1, int16(size(X,2)/2):size(X,2)) = e_xend

% Initialize Pressure
%P(1,1:size(X,2)) = const_pressure % constant pressure

% Some function to make a pressure difference along the pipe
%for i= 1:size(X,2)
    %P(1,i) = const_pressure*sin(X(i)-pi)
%end % end for

%P(1,1:int16(size(X,2)/2)) = const_pressure
%P(1,int16(size(X,2)/2):end) = const_pressure/2

P(1,1:size(X,2)) = 0 % initialization of array
P (1, 1:int16(size(X,2)/2)) = (gamma-1).*rho(1,1:int16(size(X,2)/2)).* (  e(1,1:int16(size(X,2)/2))- 0.5 .*Vx(1,1:int16(size(X,2)/2)) .*Vx(1,1:int16(size(X,2)/2))   )
P (1, int16(size(X,2)/2):end) = (gamma-1).*rho(1,int16(size(X,2)/2):end).* (  e(1,int16( size(X,2)/2 ):end)- 0.5 .*Vx(1,int16(size(X,2)/2):end) .*Vx(1,int16(size(X,2)/2):end)   )
% -------------------------------- MAIN LOOP ------------------
figure(2)
total_t = 0
for tn = 1:(size(t,2)-1) % for every time step
    dt = t(tn+1) -t(tn);
    total_t = total_t+dt
    
    %Update boundaty conditions in ghost cells
    
    %Boundary Velocity
    Vx(tn,1:(ghost_cells/2)) = Vx_x0;
    Vx(tn,(end-ghost_cells/2):end) = Vx_xend;
    
    % Boundary Density
    rho(tn,1:(ghost_cells/2)) = rho_x0;
    rho(tn,(end-ghost_cells/2):end) = rho_xend;
    
    %Boundary Pressure
    
    P(tn,1:(ghost_cells/2)) = P(1,1:(ghost_cells/2));
    P(tn,(end-ghost_cells/2):end) = P(1,(end-ghost_cells/2):end);
    
    % Boundary Specific energy 
    
    e(tn,1:(ghost_cells/2)) = e(1,1:(ghost_cells/2));
    e(tn,(end-ghost_cells/2):end) = e(1,(end-ghost_cells/2):end);
    
    % Boundary Viscosity 
    muu(tn,:) = muu(1,:); % viscosity is constant
    
    
    for xi = (ghost_cells/2+1):(size(X,2)-(ghost_cells/2)) % for every REAL position cell
        dx = X(xi+1)-X(xi);  
        dx = abs(dx);
        
        
        
        [fVx_p,fVx_m] = flux_func(tn,xi,Vx);
        S_v = -1/rho(tn,xi)*((P(tn,xi+1)-P(tn,xi-1))/(2 *dx))*dx *dt;
        Vterm = 4/3* muu(tn,xi)/rho(tn,xi)*( (Vx(tn,xi+1) -2 *Vx(tn,i) + Vx(tn,xi-1))/(dx*dx)  )*dx*dt;
        % main update equation
        %Update velocity
        Vx(tn+1,xi) = Vx(tn,xi) - dt/dx * (fVx_p - fVx_m)+ S_v;% correct solution
        
        %Update density
        [frho_p,frho_m] = flux_func(tn,xi,rho);
        rho(tn+1,xi) = rho(tn,xi) - dt/dx*(frho_p-frho_m);
        %rho(tn+1,xi) = rho(1,xi); % constant density
        
        
        %Update specific kinetic energy
        
        %Source function for energy equation
        [F_e_ip,F_e_im] = energy_flux(tn,xi,e,Vx,rho);
        S_e = - ( P(tn,xi+1)* Vx(tn,xi+1)  ) - P(tn,xi-1)*Vx(tn,xi-1)/(2*dx);
        e(tn+1,xi) = rho(tn+1,xi) * ( rho(tn,xi)*e(tn,xi) -dt/dx*(F_e_ip-F_e_im +S_e*dt)  );
        
        %Update pressure
        
        P (tn+1,xi) = (gamma-1).*rho(tn+1,xi).* (  e(tn+1,xi)- 0.5 .*Vx(tn+1,xi) .*Vx(tn+1,xi)   );
        
        
        
        
        
    end % end for every position cell
    
    scatter(X,Vx(tn,:));
    text(X(xi),Vx(tn,xi),"we are accelerated")
    timestring = sprintf("time passed = %0.5f",total_t);
    text(x_end/2,Vx_x0/2,timestring)
    pause(1e-10)
    
end % end for every time step


function [fVx_p,fVx_m] = flux_func(tn,xi,Vx)
    if (0.5*(Vx(tn,xi)+Vx(tn,xi+1)) >=0) % flow from left to right fVx_p
        fVx_p = 0.5 * (Vx(tn,xi))^2;
            
            
        else % flow from left to right
            fVx_p = 0.5 * (Vx(tn,xi+1))^2;
            
            
    end
    if 0.5*(Vx(tn,xi-1) +Vx(tn,xi)) >= 0 % if flow from left to right fVx_m
            fVx_m = 0.5 * (Vx(tn,xi-1))*(Vx(tn,xi-1) );
            
        else % if flow from right to left
            fVx_m = 0.5 * (Vx(tn,xi+1))*(Vx(tn,xi+1));
            
        end

end % end of Flux function


function [F_e_ip,F_e_im] = energy_flux(tn,xi,e,Vx,rho)
    if (0.5*(Vx(tn,xi)+Vx(tn,xi+1)) >=0) % flow from left to right fVx_p
        F_e_ip=rho(tn,xi)*Vx(tn,xi)*e(tn,xi);
        
            
            
     else % flow from left to right
         F_e_ip=rho(tn,xi+1)*Vx(tn,xi+1)*e(tn,xi+1);
            
            
            
    end
    if 0.5*(Vx(tn,xi-1) +Vx(tn,xi)) >= 0 % if flow from left to right fVx_m
        F_e_im=rho(tn,xi-1)*Vx(tn,xi-1)*e(tn,xi-1);
            
            
    else % if flow from right to left
        F_e_im=rho(tn,xi)*Vx(tn,xi)*e(tn,xi);
        
            
    end

end % end of Flux function
