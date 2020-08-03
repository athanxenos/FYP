%Cosserat Opti Test Script
%Fix tau and vary u0
clear all
close all
clc

step=-2:0.2:2;
n=length(step);
tau=[1 0 0 0];

u0=zeros(3,n);
F_err=zeros(1,n);
L_err=zeros(1,n);

u0(1,:)=step;
for i=1:n
   
    %[p] = CosseratOpti(tau(i,:),u0);
    [p] = CosseratFunc2011(tau,u0(:,i));
    
    
    hold on
    plot3(p(1,1:end-1),p(2,1:end-1),p(3,1:end-1));
    
end
    xlabel('x(m)');
    ylabel('y(m)');
    zlabel('z(m)');
    grid on
    axis([-0.25,0.25,-0.15,0.15,0,0.25]);
    title('Robot pose applying tension of 1 N in y direction for varied initial conditions')
%plot(step,F_err);
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