function [p,R,u,v] = RK2(p0,R0,u0,v0,h,tau)
%Runge Kutta 2 method of numerical integration
%Inputs:
%Initial system variables - p0,R0,u0,v0
%Tendon tension - tau
%Step size - h

%Outputs:
%System variables at next step - p,R,u,v

%Evaluate function at current step
[pd0,Rd0,ud0,vd0] = f_ode(R0,u0,v0,tau);

%Evaluate system variables at half step
R1 = R0 + Rd0*h/2;
p1 = p0 + pd0*h/2;
v1 = v0 + vd0*h/2;
u1 = u0 + ud0*h/2;
  
%Evaluate function at half step
[pd1,Rd1,ud1,vd1] = f_ode(R1,u1,v1,tau);

%Evaluate system variables at next step
R = R0 + h*Rd1;
p = p0 + h*pd1;
v = v0 + h*vd1;
u = u0 + h*ud1;
end

