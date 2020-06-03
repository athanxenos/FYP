%Cosserat Opti Test Script
clear all
close all
clc

tau=[1 0 0 0];

step=-1:0.1:1;
n=length(step);

u0=zeros(3,n);
F_err=zeros(1,n);
L_err=zeros(1,n);

u0(1,:)=step;
for i=1:n
    [p,F_err(i), L_err(i)] = CosseratOpti(tau,u0(:,i));
end

plot(step,F_err);
%Calculate arclength to check solution feasibility
%arclength = arclength(p(1,1:end-1),p(2,1:end-1),p(3,1:end-1));

% %Plot solution
% plot3(p(1,1:end-1),p(2,1:end-1),p(3,1:end-1));
% xlabel('x');
% ylabel('y');
% zlabel('z');
% grid on
% axis([-1,1,-1,1,-1,1]);
% title(['Arclength is ',num2str(arclength)])