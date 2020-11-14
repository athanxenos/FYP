function [p] = ConstantCurvature(x,y)
%ConstantCurvature
%Written by Athan Xenos

%Returns data for 2 constant curvature curves given an endpoint (x,y)
%Inputs: 
%(x,y) - Coordiante of endpoint to fit CC curves to

%Outputs: 
%p - position of CC curves

%Calculate curve midpoint
x_mid = x/2;
y_mid = y/2;

%Calculate r,theta values to describe curve
theta_mid = 2*atan(x_mid/y_mid);
r = y_mid/sin(theta_mid);

%Find coordinates of first curve
theta = linspace(0,theta_mid,10)';
x1 = r - r*cos(theta);
y1 = r*sin(theta);

%Adjuts r,theta for second curve
r2 = -r;
theta2 = linspace(theta_mid,0,10)';

%Find coordinates for second curve
x2 =   2*x_mid + r2-r2*cos(theta2);
y2 =   2*y_mid + r2*sin(theta2);

%Combine curve data
x_total = [x1;x2];
y_total = [y1;y2];

%Return coordinates for Piecewise Constant Curvature Curve
p = [x_total,y_total];
end

