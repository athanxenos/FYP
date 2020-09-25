function [s_disc] = DiscIntersect(p0,R0,v0,u0)
%Intergrates the secondary rods until they intersect the first disc
%Inputs:
%p0,R0,v0,u0 - Initial state variables of secondary rod

%Outputs:
%s_disc - Arclength where rod intersects disc

%Retrieve global variables for disc position/orientation
global d
global R_disc
global p_disc

%Set strating point before known disc location
s_disc = d(1)-0.015;

%Initialise non-zero error and counter
error_val = 1;
i=0;

tol = 0.00005;    %Tolerance in error
step = 0.0001;   %Step size (0.1mm)

%Loop until error < tolerance 
while norm(error_val) > tol
    
    %Integrate rod from 0 to s_disc
    [p,~,~,~,~] = RodODE_Eval(p0,R0,v0,u0,0, s_disc);
    
    %Check disc intersection condition and calculate error
    plane = R_disc.'*(p(end,:)-p_disc)';
    error_val = plane(end);
    
    %Increment rod arclength by step size
    s_disc = s_disc + step;
    
    %Check for infinite loop if disc intersection not found
    i=i+1;
    if i>10000
        error('Error: could not find disc intersection, check s_disc initial value/tolerance/stepsize ');
    end
end

%Subtract one step from last iteration
s_disc = s_disc - step;
end