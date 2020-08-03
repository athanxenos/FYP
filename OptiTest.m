%Cosserat Opti Test Script
%Fix u0 and vary tau
clear all
close all
clc


step=0:0.5:4;
n=length(step);
tau = zeros(n,4);
tau(:,1) = step;
u0 = [-1;0;0];
%u0=zeros(3,n);
F_err=zeros(1,n);
L_err=zeros(1,n);

%u0(1,:)=step;
for i=1:n
    if i==1
        u0=[0;0;0];
    else
        u0 = [-1;0;0];
    end
    [p] = CosseratFunc2011(tau(i,:),u0);
    %[p,F_err(i), L_err(i)] = CosseratOpti(tau,u0(:,i));
    
    
    hold on
    plot3(p(1,1:end-1),p(2,1:end-1),p(3,1:end-1));
    
end
    xlabel('x(m)');
    ylabel('y(m)');
    zlabel('z(m)');
    grid on
    axis([-0.25,0.25,-0.15,0.15,0,0.25]);
    title('Robot pose using initial conditon u0=[-1 0 0] with varied tensions')
    legend('tau = 0')
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