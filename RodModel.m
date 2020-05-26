%2D Continuum Manipulator Model
clear all
close all
clc

s=50;
d=10;
cL=55;
cR=2*s-cL;

n=length(cL);

theta = zeros(1,n);
theta_deg = zeros(1,n);
r = zeros(1,n);
L_Chord = zeros(1,n);
x = zeros(1,n);
y = zeros(1,n);



for i=1:n
    theta(i) = 2*asin((cL(i)-cR(i))/(4*d));
    theta_deg(i) = rad2deg(theta(i));
    
    r(i) = (cL(i)+cR(i))/(4*sin(theta(i)/2));
    L_Chord(i) = 2*r(i)*sin(theta(i)/2);
    
    arc_test(i) = r(i)*theta(i);
    r_test(i) = s/theta(i);
    theta_test(i) = s/r(i);
    
    x(i) = r(i)*(1-cos(theta(i)));
    y(i) = r(i)*sin(theta(i));
end

hold on

for i=1:n
    [xArc, yArc] = plotArc(x(i),r(i));
    plot(xArc,yArc);
    axis([0,50,0,50]);
end