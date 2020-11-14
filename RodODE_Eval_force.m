function [p,R,n,m,s] = RodODE_Eval_force(p0,R0,n0,m0,L_start,L_finish)
%Function that evaluates ODE's from L_start to L_finish
%Written by Athan Xenos

%Inputs:
% p0,R0,n0,m0 - Initial state variables of rod 
% L_start - starting length
% L_finish - final length

%Outputs:
% p,R,n,m - State variables along rod length
% s - Vector of 100 elements from L_start to L_finish

%Create evenly spaced vector of length 100 
L = linspace(L_start,L_finish);

%Setup initial values for ODE
y0 = [p0 ; reshape(R0,9,1); n0; m0];

%Solve ODE from L_start to L_finish
[s,y] = ode45(@rod_ode_force, L, y0);

%Extract solution curve values
n = length(s);

%Define curve position
p= [y(:,1) y(:,2) y(:,3)];

%Extract state variables along rod
R = reshape(y(:,4:12)',3,3,n);
n = [y(:,13) y(:,14) y(:,15)];
m = [y(:,16) y(:,17) y(:,18)];
end