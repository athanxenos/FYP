%Cosserat Opti Test Script
%Runs Cosserat Func Multiple Times to optimise initial guess
%Fix tau and vary u0
clear all
close all
clc

%Set tau input
tau=[2.5 0 0 0];

%Define u0 step size
step=-2:0.1:2;
n=length(step);

%Initialise matrices/vector
u0=zeros(3,n);
F_err=zeros(1,n);
L_err=zeros(1,n);

%Set u0 variable
u0(1,:)=step;

%Loop through each initial condition
for i=1:n
   
    %Evaluate Cosserat Function for each initial guess
    [p,F_err(i),L_err(i)] = CosseratFunc2011(tau,u0(:,i));
    
    %Plot each rod position
    hold on
    plot3(p(1,1:end-1),p(2,1:end-1),p(3,1:end-1));
    
end

%     xlabel('x(m)');
%     ylabel('y(m)');
%     zlabel('z(m)');
%     grid on
%     axis([-0.25,0.25,-0.15,0.15,0,0.25]);
%     title('Robot pose applying tension of 1 N in y direction for varied initial conditions')
    
figure    
plot(step,F_err);

xlabel('ux');
ylabel('F_error');
title('F error for fixed tau');
grid on

figure
plot(step,L_err);

xlabel('ux');
ylabel('L_error');
title('L error for fixed tau');
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