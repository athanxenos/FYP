%2D Continuum Manipulator Model
clear all
close all
clc

s=50;
d=10;
qL=50.1:1:s+d*pi;
qR=2*s-qL;

n=length(qL);

theta = zeros(1,n);
theta_deg = zeros(1,n);
r = zeros(1,n);

x = zeros(1,n);
y = zeros(1,n);

for i=1:n
    theta(i) = (qL(i)-s)/d;
    theta_deg(i) = rad2deg(theta(i));
    
    r(i) = s/theta(i);
    
    x(i) = r(i)*(1-cos(theta(i)));
    y(i) = r(i)*sin(theta(i));
end

hold on

for i=1:n
    [xArc, yArc] = plotArc(x(i),r(i));
    plot(xArc,yArc);
    axis([0,50,0,50]);
end

sTest = 50;
rTest = linspace(s/pi,10*s);
thetaTest = sTest./rTest;
xTest = rTest.*(1-cos(thetaTest));
yTest = rTest.*sin(thetaTest);

%plot(xTest,yTest);

for i=1:length(rTest)
    [xArc, yArc] = plotArc(xTest(i),rTest(i));
    %plot(xArc,yArc);
    axis([0,50,0,50]);
end
