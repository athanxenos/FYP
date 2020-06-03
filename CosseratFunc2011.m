%Cosserat Model
clear all
close all
clc 

%Rod Parameters
L = 1; %Arclength of rod (m)
ds = 0.01; %Step size
s = 0:ds:L; %arclength parameter (m)
n = length(s);


%Initial Conditions
R0 = eye(3); %Initial rod orientation at base
p0 = [0;0;0]; %Inital rod position at base
v0=[0;0;1]; %Guess initial conditions for v,u
u0=[-1;0;0];


%Setup initial iteration
R=zeros(3,3,n);
p=zeros(3,n);
v=zeros(3,n);
u=zeros(3,n);
n_rod=zeros(3,n);
m=zeros(3,n);


%Assign initial values
R(:,:,1)=R0;
p(:,1)=p0;
v(:,1)=v0;
u(:,1)=u0;



for j=1:n-1
    
    %n_rod(:,j) = R(:,:,j)*K_se*(v(:,j)-v_ref);
    %m(:,j) = R(:,:,j)*K_bt*(u(:,j)-u_ref);
%     
%     [pd,Rd,ud,vd] = f_ode(p(:,j),R(:,:,j),u(:,j),v(:,j));
%     
%     if j<101
%      R(:,:,j+1) = R(:,:,j) + Rd*ds;
%      p(:,j+1) = p(:,j) + pd*ds;
%      v(:,j+1) = v(:,j) + vd*ds;
%      u(:,j+1) = u(:,j) + ud*ds;
%     end
    
[p(:,j+1),R(:,:,j+1),u(:,j+1),v(:,j+1)] = RK2(p(:,j),R(:,:,j),u(:,j),v(:,j),ds); 

end



arclength(p(1,:),p(2,:),p(3,:))

plot3(p(1,:),p(2,:),p(3,:));
xlabel('x');
ylabel('y');
zlabel('z');
grid on
axis([-1,1,-1,1,-1,1]);