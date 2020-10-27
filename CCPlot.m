%Results Plot
clear variables
clear global
close all
clc 

data = load('M_0.2.mat');
pb = data.pb;

%Plot solution for central backbone
hold on
plot(pb(:,2),pb(:,3),'b');


[p] = ConstantCurvature(pb(end,2),pb(end,3));
plot(p(:,1),p(:,2),'b.');