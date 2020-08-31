function [p,R,v,u] = RodIntegrate(p0,R0,v0,u0,L)
%Function that evaluates ODE's for input guess and returns residual when
%compared to known boundary conditions

%Inputs:
% Initial state of rod described by - p0,R0,v0,u0
% Length of rod to integrate - L
%Outputs:
% State variables along rod length - p,R,v,u

%Setup initial values for ODE
y0 = [p0 ; reshape(R0,9,1); v0; u0];

%Solve ODE from 0 to L
[s,y] = ode45(@f_secondary, [0 L], y0);

%Extract solution curve values
n = length(s);
px = y(:,1);
py = y(:,2);
pz = y(:,3);

%Define curve as global
p= [px py pz];

%Extract state variables at endpoint s=L
R = reshape(y(:,4:12)',3,3,n);
v = [y(:,13) y(:,14) y(:,15)];
u = [y(:,16) y(:,17) y(:,18)];
end