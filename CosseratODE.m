function [yd] = CosseratODE(s,y)
%CosseratODE
%Written by Athan Xenos

%Function to evaluate ODE system at one step for Cosserat Rod model
%Inputs:
%y (p,R,n,m) - Systems variables described as state vector  
%s - Rod length, not used in equations, used for ode45 input 

%Outputs:
%yd (pd,Rd,nd,md) - State vector derivative

%Global variables determined by robot structure
global K_se
global K_bt
global v_ref

%Unpack state vector
R = reshape(y(4:12),3,3);
n = y(13:15);
m = y(16:18);

%Derivative calculation
v = K_se^-1*R'*n + v_ref;
u = K_bt^-1*R'*m;
f = [0;0;0];

pd = R*v;
Rd = R*hat(u);
nd = -f;
md = -cross(pd,n);

%Repack derivative vector
yd = [pd ; reshape(Rd,9,1); nd; md];
end