%Cosserat Opti Test Script
%Runs Cosserat Func Multiple Times to optimise initial guess
%Fix u0 and vary tau
clear all
close all
clc

%Set initial guess input for u
u0 = [-0.5;0;0];

%Define tau step size
step=0:0.5:5;
n=length(step);
tau = zeros(n,4);
tau(:,1) = step;

%Initialise matrices/vector
F_err=zeros(1,n);
L_err=zeros(1,n);


for i=1:n
%     if i==1
%         u0=[0;0;0];
%     else
%         u0 = [-1;0;0];
%     end
    [p,F_err(i), L_err(i)] = CosseratFunc2011(tau(i,:),u0);
    
    
    hold on
    plot3(p(1,:),p(2,:),p(3,:));
    
end
    xlabel('x(m)');
    ylabel('y(m)');
    zlabel('z(m)');
    grid on
    axis([-0.25,0.25,-0.15,0.15,0,0.25]);
    title('Robot pose using initial conditon u0=[-1 0 0] with varied tensions')
    

    
figure    
plot(step,F_err);

xlabel('tau');
ylabel('F_error');
title('F error for fixed u0');
grid on

figure
plot(step,L_err);

xlabel('tau');
ylabel('L_error');
title('L error for fixed u0');
grid on


%Calculate arclength to check solution feasibility
arclength = arclength(p(1,1:end-1),p(2,1:end-1),p(3,1:end-1));

% %Plot solution
% plot3(p(1,1:end-1),p(2,1:end-1),p(3,1:end-1));
% xlabel('x');
% ylabel('y');
% zlabel('z');
% grid on
% axis([-1,1,-1,1,-1,1]);
% title(['Arclength is ',num2str(arclength)])