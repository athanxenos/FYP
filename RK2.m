function [p,R,u,v,n,m] = RK2(p0,R0,u0,v0,tau)
%Runge Kutta 2: Method of numerical integration

%Inputs:
%Initial system variables - p0,R0,u0,v0
%Tension vector - tau

%Outputs:
%System variables at next step (ds along robot) - p,R,u,v

%Import global variables
global K_se
global K_bt
global v_ref
global u_ref
global ds  %Step size

%Evaluate function at current step
[~,Rd0,ud0,vd0] = f_ode(R0,u0,v0,tau);

%Evaluate system variables at half step (ds/2)
R1 = R0 + Rd0*ds/2;
v1 = v0 + vd0*ds/2;
u1 = u0 + ud0*ds/2;
  
%Evaluate function at half step (ds/2)
[pd1,Rd1,ud1,vd1] = f_ode(R1,u1,v1,tau);

%Evaluate system variables at next step (ds)
R = R0 + ds*Rd1;
p = p0 + ds*pd1;
v = v0 + ds*vd1;
u = u0 + ds*ud1;

%Evaluate n,m at next step
n = R*K_se*(v-v_ref);
m = R*K_bt*(u-u_ref);
end