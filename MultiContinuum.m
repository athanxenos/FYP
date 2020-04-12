%Multisection Continuum Robot Modelling
%Jones Kinematic Model
%Section 1
clear all
close all
clc

%Actuator Parameters
L1=10;
L2=10;
L3=15;

n=3; %Number of cable guide per section
d=2; %Radius of trunk

g = L1^2+L2^2+L3^2-L1*L2-L2*L3-L1*L3;

%Arc Length
%s = (n*d*(L1+L2+L3)/sqrt(g))*asin(sqrt(g)/(3*n*d));
s=15;

%Curvature
%k = 2*sqrt(g)/(d*(L1+L2+L3));
k=-0.1;

%Direction of Curvature Angle (radians)
%phi = atan((sqrt(3)/3) *((L3+L2-2*L1)/(L2-L3)));
%phi_deg = rad2deg(phi);
phi = pi/3;

s_step = 0:1:s;
s_step(end+1) = s;

for i=1:length(s_step)
    x(i) = sin(phi)*(1-cos(s_step(i)*k))/k;
    y(i) = -cos(phi)*(1-cos(s_step(i)*k))/k;
    z(i) = sin(s_step(i)*k)/k;
end

T=TMatrix(s,k,phi);

plot3(x,y,z);
grid on
%axis equal