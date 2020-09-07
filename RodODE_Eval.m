function [p,R,v,u,s] = RodODE_Eval(p0,R0,v0,u0,L_start,L_finish)
%Function that evaluates ODE's from L_start to L_finish

%Inputs:
% p0,R0,v0,u0 - Initial state variables of rod 
% L_start - starting length
% L_finish - final length

%Outputs:
% p,R,v,u - State variables along rod length
% s - Vector of 100 elements from L_start to L_finish

%Create evenly spaced vector of length 100 
L = linspace(L_start,L_finish);

%Setup initial values for ODE
y0 = [p0 ; reshape(R0,9,1); v0; u0];

%Solve ODE from L_start to L_finish
[s,y] = ode45(@rod_ode, L, y0);

%Extract solution curve values
n = length(s);

%Define curve position
p= [y(:,1) y(:,2) y(:,3)];

%Extract state variables along rod
R = reshape(y(:,4:12)',3,3,n);
v = [y(:,13) y(:,14) y(:,15)];
u = [y(:,16) y(:,17) y(:,18)];
end