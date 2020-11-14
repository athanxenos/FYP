%CCPlot
%Written by Athan Xenos

%Script that fits Constant Curvature Curves to simulation data from MultiBackboneOpti
%Plots resulting comparison
clear variables
clear global
close all
clc 

%Load simulation data
data2 = load('M_0.2.mat');
pb2 = data2.pb;
data4 = load('M_0.4.mat');
pb4 = data4.pb;
data6 = load('M_0.6.mat');
pb6 = data6.pb;
data8 = load('M_0.8.mat');
pb8 = data8.pb;

%Plot solution for central backbone
hold on

%Fit Constant Curvature curves to data
p2 = ConstantCurvature(pb2(end,2),pb2(end,3));
p4 = ConstantCurvature(pb4(end,2),pb4(end,3));
p6 = ConstantCurvature(pb6(end,2),pb6(end,3));
p8 = ConstantCurvature(pb8(end,2),pb8(end,3));

%Plot data and CC curves
plot(pb2(:,2),pb2(:,3),'b',pb4(:,2),pb4(:,3),'r',pb6(:,2),pb6(:,3),'m',pb8(:,2),pb8(:,3),'k');
plot(p2(:,1),p2(:,2),'bo',p4(:,1),p4(:,2),'ro',p6(:,1),p6(:,2),'mo',p8(:,1),p8(:,2),'ko');

%Calculate distance and error between curves
[k2,dist2] = dsearchn(pb2(:,2:3),p2);
error2 = norm(dist2);
[k4,dist4] = dsearchn(pb4(:,2:3),p4);
error4 = norm(dist4);
[k6,dist6] = dsearchn(pb6(:,2:3),p6);
error6 = norm(dist6);
[k8,dist8] = dsearchn(pb8(:,2:3),p8);
error8 = norm(dist8);

avg2 = mean(dist2(2:end-1));
avg4 = mean(dist4(2:end-1));
avg6 = mean(dist6(2:end-1));
avg8 = mean(dist8(2:end-1));

%Graph labels
xlabel('y (m)');
ylabel('z (m)');
grid on
axis([0,0.3,0,0.3]);
title('Cosserat Rod Model (2 Point Moments) vs Constant Curvature Model');
legend('Cosserat Pose 1','Cosserat Pose 2','Cosserat Pose 3','Cosserat Pose 4','Constant Curvature Arcs');