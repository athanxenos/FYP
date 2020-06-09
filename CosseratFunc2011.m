%Cosserat Model
clear all
close all
clc 

global K_se
global K_bt
global r
global v_ref
global u_ref
global ud_ref
global vd_ref

%Rod Parameters
L = 0.25; %Arclength of rod (m)
ds = 0.01; %Step size
s = 0:ds:L; %arclength parameter (m)
n = length(s);

rad = 0.0005; %Radius of central rod (m)(approx 0.5mm)

%Material parameters
Area = pi*rad^2; %Area of cross section (m^2)
EY = 210*10^9; %Young's Modulus (Pa)
GY = 80*10^9;  %Shear Modulus (Pa) G=E/2(1+v)

%Second moments of area of cross section
Ixx = pi*rad^4/4;
Iyy = pi*rad^4/4;
Izz = Ixx + Iyy; %Polar moment of inertia

%Stiffness Matrices
K_se = diag([GY*Area GY*Area EY*Area]); %Shear/extension stiffness matrix
K_bt = diag([EY*Ixx EY*Iyy EY*Izz]); %Bending/torsion stiffness matrix

%Tendon Parameters
n_t = 4; %Number of tendons
r_t = 0.01; %Radial location of tendons (m)(approx 10mm)

%x,y locations of tendons in rod cross section
x_t = [0 -r_t 0 r_t];
y_t = [r_t 0 -r_t 0];
r = [x_t; y_t; 0 0 0 0]; %position vectors of tendons

%Reference Parameters
%Linear/angular rate of change of frame in reference state
v_ref = [0;0;1];
u_ref = [0;0;0];
ud_ref = [0;0;0];
vd_ref = [0;0;0];

%Initial Conditions
R0 = eye(3); %Initial rod orientation at base
p0 = [0;0;0]; %Inital rod position at base

%Guess initial conditions for v,u
v0=[0;0;1];
u0=[0;0;0];

%Tension Input
tau = [0 0 0 0]; %Tension for each tendon

%Setup initial iteration
R=zeros(3,3,n);
p=zeros(3,n);
v=zeros(3,n);
u=zeros(3,n);


n_rod=zeros(3,n);
m=zeros(3,n);

pid = zeros(3,n_t);
F_tendon = zeros(3,n_t);
L_tendon = zeros(3,n_t);

%Assign initial values
R(:,:,1)=R0;
p(:,1)=p0;
v(:,1)=v0;
u(:,1)=u0;

n_rod(:,1) = R0*K_se*(v0-v_ref);
m(:,1) = R0*K_bt*(u0-u_ref);

%Iterate solution from s=0 to s=L using RK2
for j=1:n
    
 
[p(:,j+1),R(:,:,j+1),u(:,j+1),v(:,j+1)] = RK2(p(:,j),R(:,:,j),u(:,j),v(:,j),ds,tau); 

n_rod(:,j+1) = R(:,:,j+1)*K_se*(v(:,j+1)-v_ref);
m(:,j+1) = R(:,:,j+1)*K_bt*(u(:,j+1)-u_ref);

end



for i=1:n_t
    pid(:,i) = R(:,:,n)*(hat(u(:,n))*r(:,i)+v(:,n));
    F_tendon(:,i) = -tau(i)*pid(:,i)/norm(pid(:,i));
    L_tendon(:,i) = -tau(i)*hat(R(:,:,n)*r(:,i))*pid(:,i)/norm(pid(:,i));
end

F_sum = sum(F_tendon,2);
L_sum = sum(L_tendon,2);

F_error = norm(F_sum-n_rod(:,n));
L_error = norm(L_sum-m(:,n));

%F_err_abs = norm(F_sum-(n_rod(:,n-1)-n_rod(:,n+1)))
%L_err_abs = norm(L_sum-(m(:,n-1)-m(:,n+1)))


%Calculate arclength to check solution feasibility
arclength = arclength(p(1,1:end-1),p(2,1:end-1),p(3,1:end-1))

%Plot solution
plot3(p(1,1:end-1),p(2,1:end-1),p(3,1:end-1));
xlabel('x');
ylabel('y');
zlabel('z');
grid on
axis([-0.25,0.25,-0.25,0.25,-0.25,0.25]);
title(['Arclength is ',num2str(arclength)])