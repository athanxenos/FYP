function [p,R,u,v] = RK2(p0,R0,u0,v0,h)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
[pd0,Rd0,ud0,vd0] = f_ode(p0,R0,u0,v0);

     R1 = R0 + Rd0*h/2;
     p1 = p0 + pd0*h/2;
     v1 = v0 + vd0*h/2;
     u1 = u0 + ud0*h/2;
     
[pd1,Rd1,ud1,vd1] = f_ode(p1,R1,u1,v1);

     R = R0 + h*Rd1;
     p = p0 + h*pd1;
     v = v0 + h*vd1;
     u = u0 + h*ud1;
end

