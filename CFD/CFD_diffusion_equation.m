clc
clear
%ggg T - temp ;  t- time
% T[xi](n)   - temperature in cell with index i, at time cell t = n

% ∂T(x,t)/ ∂t  = (T[xi](n+1) -T[xi](n) )/ dt    numerical derivative 
% ∂∂T(x,t)/ ∂x"2 = T[xi+1](n) -2 T[xi](n) + T[i-1](n)/ dx"2 

% where : ∂∂ - second order partial derivative
% d - delta x between cells
% T(x,t) - continuous function of temperature with respect to possition x
% and time t
% ∂∂T(x,t)/ ∂x"2   - second order derivative with respect to position x
% T[xi](n)  - temperature value at descreete cell with position [xi] and
% time cell t=n. 

% general continious equation :   ∂T(x,t)/ ∂t  = α* (∂∂T(x,t)/ ∂x"2)
% where α - ?

% discreete equation : (T[xi](n+1) -T[xi](n) )/ dt  = (T[xi+1](n) -2 T[xi](n) + T[i-1](n) )/ dx"2
% where dt= t(n+1) - t(n)  difference in time between time cells n+1 and n

% from discreete equation we derive :
%T(xi)[n+1] = T(xi)[n] +α*dt *  ( T[xi+1](n) - 2 T[xi](n) + T[xi-1](n)  )/dx"2


% DEFINE parameters
x0 =0
xend = 0.10
alpha = 19e-5; %thermal conductivity
dt = 2; %timestep
tMax = 10000
real_cells = 10
%Define T at initial time step ------------------
%boundary condition temperatures
T_x0 = 353
T_xend = 293

T(1,1:2) = T_x0
T(1,3:12) = T_xend 

%COLUMNS - position
% ROWS - timestep


%DEFINE X position cells ------------------
%dif = (xend-x0)/real_cells
%X = linspace(x0-dif,xend+dif,real_cells+2)

dX = (xend -x0)/(real_cells-2)

for i=1:real_cells+2
    X(i) = x0 +(i-3/2)*dX
end

% Timestep vector -----------------------------
time_cells_number = int16( tMax / dt)
t = linspace(0,tMax,time_cells_number);



% Main loop

for tn =1:time_cells_number-1 % for every timedtep
    dt = t(tn+1) -t(tn)
    %Now set boundaty conditions in ghost cells
    T(tn,1) = T_x0;
    T(tn,end) = T_xend;
    for xi = 2:(size(X,2)-1) % for every REAL position cell
        dx = X(xi+1)-X(xi);  
        %%T(xi)[n+1] = T(xi)[n] +α*dt *  ( T[xi+1](n) - 2 T[xi](n) + T[xi-1](n)  )/dx"2
        % update Temperature 
        T(tn+1,xi) = T(tn,xi) +alpha *dt * ( T(tn,xi+1) - 2* T(tn,xi) + T(tn,xi-1) )/dx*dx;
    end
    
   
    
    
    
end