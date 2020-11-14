%Results Plot
clear variables
clear global
close all
clc 

%Import experimental data and fit 2 constant curvature arcs to endpoint
x20 = [0.11,2.47,9.5,18.88,30.39,40.88,47.86,52.13,53.85]/1000;
y20 = [0,35.7,71.39,103.97,138.58,174.75,208.95,243.69,278.51]/1000;
p20 = ConstantCurvature(x20(end),y20(end));

x40 = [0.14,4.57,18.11,35.7,57.13,77.12,90.21,98.19,101.5]/1000;
y40 = [0,35.82,69.71,98.23,127.71,159.81,192.08,226.3,261.01]/1000;
p40 = ConstantCurvature(x40(end),y40(end));

x60 = [0.18,6.76,25.93,50.18,79.78,107.35,125.5,136.91,141.92]/1000;
y60 = [0,35.74,66.43,89.06,110.89,136.78,166.29,199.55,234.13]/1000;
p60 = ConstantCurvature(x60(end),y60(end));

x80 = [0.29,9.23,34.34,63.77,98.65,131.17,153.28,167.37,174.76]/1000;
y80 = [0,35.07,61.17,75.82,87.13,105.39,132.63,165.03,198.92]/1000;
p80 = ConstantCurvature(x80(end),y80(end));

[k2,dist2] = dsearchn(p20,[x20',y20']);
error2 = norm(dist2);
[k4,dist4] = dsearchn(p40,[x40',y40']);
error4 = norm(dist4);
[k6,dist6] = dsearchn(p60,[x60',y60']);
error6 = norm(dist6);
[k8,dist8] = dsearchn(p80,[x80',y80']);
error8 = norm(dist8);

avg2 = mean(dist2(1:end-1));
avg4 = mean(dist4(1:end-1));
avg6 = mean(dist6(1:end-1));
avg8 = mean(dist8(1:end-1));

%Plot experiment data vs CC curves
hold on
plot(x20,y20,'bd',x40,y40,'rd',x60,y60,'md',x80,y80,'kd');
plot(p20(:,1),p20(:,2),'b',p40(:,1),p40(:,2),'r',p60(:,1),p60(:,2),'m',p80(:,1),p80(:,2),'k');

%Graph labels
xlabel('y (m)');
ylabel('z (m)');
grid on
axis([0,0.3,0,0.3]);

title('Constant Curvature Model vs Experimental Data')
legend('Experiment 1 Data','Experiment 2 Data','Experiment 3 Data','Experiment 4 Data','Two CC Arcs','Two CC Arcs','Two CC Arcs','Two CC Arcs','Location','SouthEast');
