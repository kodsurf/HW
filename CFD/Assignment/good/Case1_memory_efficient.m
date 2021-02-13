%CFD ASSIGNMENT
clc
clear
% CFD 

% CASE 1

%ASSUMPTIONS: 
% 1) CONSTANT VELOCITY THROUGH OUT THE PIPE
% 2) CONSTANT THERMAL PROPERTIES
% 3) INCOPMRESSBLE WATER FLOW
% 4) ADIABATIC WALLS 
% ------------ INITIAL PARAMETERS ------------------
PSEUDO_REAL_TIME = 0 % live plot  1 - on  0 - off
MEMORY_EFFICIENT = 1 % MORE FAST if 1 / Save every time step temperature disctibution 0
USE_CUSTOM_DELTA_T = 0 % if 1 - user can specify dt manually  0 - calcualted according to CFL conditions 
% Number of cells
real_cells = 150
ghost_cells = 2 % should be even

% ------ Pipe length 
x0 = 0
x_end =0.04
% -- Pipe area
A = 0.002 %mm
d = 0.002

% -- Fluid velocity Vx
U_constant = 0.004 % m/s


% --Density 
const_rho = 1000; % kg/m"3
%rho_x0 = 1
%rho_xend = 0.125

% Thermal conductivity
k =0.6 % W/K
% Specific heat capacity
c = 4200 % J (kg*K)
%Convective heat transfer coefficient
h = 3000 % W/(m"2*K)

% Time of computation and timestep 
user_defined_dt = 1e-5*100
tMax = 5



% Boundary conditions for TEMPERATURE (watter)
T_x0 = 400 %K
T_xend = 293%K

T_air = 293  %270 %K
T_wall = 400% K






% ---------------------- INITIALIZATION

% Initialize X cells -------------------
dX = (abs(x0 -x_end))/(real_cells);

for i=1:real_cells+ghost_cells
    %X(i) = x0 +(i-3/2)*dX;
    X(i) = x0+(i-(ghost_cells/2+0.5))*dX;
end

% Initialize Time cells

if USE_CUSTOM_DELTA_T==1
    dt = user_defined_dt
else
    %*k*dt/(dx) <1
    dt = dX/(k) % at least 
    dt = 0.8*dt % multiply by factor so that dt < min required dt
end
time_cells_number = int16( tMax / dt)
time_cells_number = double(time_cells_number)
%t = linspace(0,tMax,time_cells_number);
t(1) = 0;
for i=2:time_cells_number
    
    t(i) = dt*i;
end

% Initialize TEMPERATURE

Tx(1,1:(ghost_cells/2)) = T_x0
Tx(1,2:(size(X,2))) = T_xend;

Tx_old(1:ghost_cells/2) = T_x0
Tx_old(2:(size(X,2))) = T_xend;

Tx_new(1:ghost_cells/2) = T_x0
Tx_new(2:(size(X,2))) = T_xend;

% Initialize Velocity

U(1,1:int16(size(X,2))) = U_constant



% -------------------------------- MAIN LOOP ------------------
figure(2)
%hold on
total_t = 0
tn=1
%for tn = 1:(size(t,2)-1) % for every time step
while total_t <= tMax
    
    if MEMORY_EFFICIENT ==1
        
    else
        dt = t(tn+1) -t(tn);
        tn = tn +1;
    end
    total_t = total_t+dt
    
    %Update boundaty conditions in ghost cells
    
    %IDENTIFY BOUNDARY CELL
    flow_travel_distance = total_t*U_constant;
    
    %identify index of boundary cell
    %boundary_cell_index = (flow_travel_distance/dX)
    boundary_cell_index =int16(0.5+flow_travel_distance/dX);
    if boundary_cell_index==0
        boundary_cell_index = 1;
    end
    boundary_cell_index = boundary_cell_index % debug
        
    
    %Boundary velocity
    % CONSTANT 
    %U(tn,:) = U_constant;
    
    if MEMORY_EFFICIENT == 1
        %Boundary Temperature ( AT AIR-WATER BORDER)
        Tx_old(1:(ghost_cells/2)) = T_x0; % INFLOW TEMPERATURE
        %Tx(tn,boundary_cell_index:end) = T_xend;
        Tx_old(boundary_cell_index) = T_xend;
        Tx_old(boundary_cell_index+1:end) = T_air; % temperature of air
    else
        Tx(tn,1:(ghost_cells/2)) = T_x0; % INFLOW TEMPERATURE
        %Tx(tn,boundary_cell_index:end) = T_xend;
        Tx(tn,boundary_cell_index) = T_xend;
        Tx(tn,boundary_cell_index+1:end) = T_air; % temperature of air
    end
    


    
    for xi = (ghost_cells/2+1):boundary_cell_index % for every REAL position TO WHICH FLOW REACHED
        dx = X(xi+1)-X(xi);  
        dx = abs(dx);
        
        
        

        
        %Update Temperature
        if xi>=boundary_cell_index;
            %If boundary cell - dont evaluate diffusive term (because it
            %uses T(x+1). There is no heat transwer benween air/water. Thus
            %using T(x+1) would be air cell. Thus this would have been
            %wrong
            if MEMORY_EFFICIENT == 1 %Save memory
                Tx_new(xi) = 0*k*dt/(dx*dx*const_rho*c)*(Tx_old(xi+1) -2*Tx_old(xi) + Tx_old(xi-1)) +Tx_old(xi) - dt*U_constant*(Tx_old(xi)-Tx_old(xi-1))/dx;%+4*h*dt*(T_wall-Tx_old(xi))/(c*const_rho*d);
                
            else % Save all progress
                %Tx(tn+1,xi) = 0*k*dt/(dx*dx*const_rho*c)*(Tx(tn,xi+1) -2*Tx(tn,xi) + Tx(tn,xi-1)) +Tx(tn,xi) - dt*U(tn,xi)*(Tx(tn,xi)-Tx(tn,xi-1))/dx+4*h*dt*(T_wall-Tx(tn,xi))/(c*const_rho*d);
                Tx(tn+1,xi) = 0*k*dt/(dx*dx*const_rho*c)*(Tx(tn,xi+1) -2*Tx(tn,xi) + Tx(tn,xi-1)) +Tx(tn,xi) - dt*U_constant*(Tx(tn,xi)-Tx(tn,xi-1))/dx;%+4*h*dt*(T_wall-Tx(tn,xi))/(c*const_rho*d);
            end
            
            
        else
            if MEMORY_EFFICIENT == 1 %Save memory
                Tx_new(xi) = k*dt/(dx*dx*const_rho*c)*(Tx_old(xi+1) -2*Tx_old(xi) + Tx_old(xi-1)) +Tx_old(xi) - dt*U_constant*(Tx_old(xi)-Tx_old(xi-1))/dx;%+4*h*dt*(T_wall-Tx_old(xi))/(c*const_rho*d);
            else % Save progress
                %Tx(tn+1,xi) = k*dt/(dx*dx*const_rho*c)*(Tx(tn,xi+1) -2*Tx(tn,xi) + Tx(tn,xi-1)) +Tx(tn,xi) - dt*U(tn,xi)*(Tx(tn,xi)-Tx(tn,xi-1))/dx+4*h*dt*(T_wall-Tx(tn,xi))/(c*const_rho*d);
                Tx(tn+1,xi) = k*dt/(dx*dx*const_rho*c)*(Tx(tn,xi+1) -2*Tx(tn,xi) + Tx(tn,xi-1)) +Tx(tn,xi) - dt*U_constant*(Tx(tn,xi)-Tx(tn,xi-1))/dx;%+4*h*dt*(T_wall-Tx(tn,xi))/(c*const_rho*d);
            end
        end
        
        
        
        
    end % end for every position cell
    if PSEUDO_REAL_TIME == 1
        if MEMORY_EFFICIENT == 1
            scatter(X,Tx_new);
        else
            scatter(X,Tx(tn,:));
        end
    title('Case 2: Adiabatic walls')
    xlabel('X (m)') 
    ylabel('T (K)')
    timestring = sprintf("time passed = %0.2f s",round(total_t,2));
    text(0,300,timestring)
    pause(1e-10)
    end

    
    if MEMORY_EFFICIENT == 1
        Tx_old = Tx_new
    end
    
end % end for every time step

figure(3)

if MEMORY_EFFICIENT == 1
    scatter(X,Tx_new);
else
    scatter(X,Tx(tn,:));
end
title('Case 1: Adiabatic walls')
xlabel('X (m)') 
ylabel('T (K)')
velocity_string = sprintf("u= %0.5f m/s",U_constant);
timestring = sprintf("Time passed = %0.2f s",round(total_t,2));
text(0,300,timestring)
text(0,310,velocity_string)
