function [yd] = f_secondary(s,y)
%Function to evaluate ODE system at one step for secondary backbones

%Inputs:
%Systems variables described as state vector - y (p,R,v,u)
%Rod length, not used in equations, used for ode45 input - s

%Outputs:
%State vector derivative - yd (pd,Rd,vd,ud)

%Global variables determined by robot structure
global K_se
global K_bt
global v_ref


%Unpack state vector
R = reshape(y(4:12),3,3);
v = y(13:15);
u = y(16:18);

%Derivative calculation
pd = R*v;
Rd = R*hat(u);

vd = -K_se^-1*(hat(u)*K_se*(v-v_ref));
ud = -K_bt^-1*(hat(u)*K_bt*u+hat(v)*K_se*(v-v_ref));

yd = [pd ; reshape(Rd,9,1); vd; ud];
end