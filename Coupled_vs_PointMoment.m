%Results Plot
clear variables
clear global
close all
clc 

data = load('Tendon_M0.2.mat');
pb_tendon = data.pb;

data2 = load('PointMoment_M0.2.mat');
pb_moment = data2.pb;

plot(pb_tendon(:,2),pb_tendon(:,3),'b',pb_moment(:,2),pb_moment(:,3),'r');


%Graph labels
xlabel('y (m)');
ylabel('z (m)');
grid on
axis([0,0.1,0,0.3]);

title('Coupled Tendon Model vs Point Moment Model for 0.2Nm Moment');
legend('Coupled Tendon Model','Point Moment Model');