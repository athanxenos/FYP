%Cosserat Model Script based on 2011 paper
%Uses RK2 and  f_ode functions to solve model
clear all
close all
clc 

%Define global variables for model
global K_se
global K_bt
global r
global v_ref
global u_ref
global ud_ref
global vd_ref
global ds

%Rod Parameters
L = 0.25; %Arclength of rod (m)
rad = 0.0005; %Radius of central rod (m)(approx 0.5mm)
Area = pi*rad^2; %Area of cross section (m^2)

%Discretisation Parameters
ds = 0.01; %Step size
s = 0:ds:L; %arclength parameter (m)
n = length(s);

%Material parameters
EY = 210*10^9; %Young's Modulus (Pa)
GY = 80*10^9;  %Shear Modulus (Pa) G=E/2(1+v) (v=0.3125)

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
r = [x_t; y_t; 0 0 0 0]; %Position vectors of tendons in body frame

%Reference Parameters
%Linear/angular rate of change of frame in reference state
v_ref = [0;0;1];
u_ref = [0;0;0];
ud_ref = [0;0;0];
vd_ref = [0;0;0];

%/////////// Model Variables ////////////
%Guess initial conditions for v,u
v0=[0;0;1]; %Linear rate of change of frame
u0=[-1;0;0]; %Angular rate of change of frame

%Tension Input
tau = [1 0 0 0]; %Tension for each tendon (N)
%////////////////////////////////////////

%Initial Conditions
R0 = eye(3); %Initial rod orientation at base
p0 = [0;0;0]; %Inital rod position at base
n0= R0*K_se*(v0-v_ref); %Initial internal force
m0= R0*K_bt*(u0-u_ref); %Initial internal moment

%Initialise matrices/vectors for system parameters
R=zeros(3,3,n);
p=zeros(3,n);
v=zeros(3,n);
u=zeros(3,n);

n_rod=zeros(3,n);
m=zeros(3,n);

%Assign initial values at base
R(:,:,1)= R0;
p(:,1)= p0;
v(:,1)= v0;
u(:,1)= u0;
n_rod(:,1) = n0;
m(:,1) = m0;

%Initialise key variables
pid = zeros(3,n_t);
F_tendon = zeros(3,n_t);
L_tendon = zeros(3,n_t);


%Iterate solution from s=0 to s=L using RK2

for j=1:n-1  
%Pass system variables into RK2 function to update for each step
[p(:,j+1),R(:,:,j+1),u(:,j+1),v(:,j+1),n_rod(:,j+1),m(:,j+1)] = RK2(p(:,j),R(:,:,j),u(:,j),v(:,j),tau);
end


%Iterate through each tendon
for i=1:n_t
    pid(:,i) = R(:,:,n)*(hat(u(:,n))*r(:,i)+v(:,n));
    F_tendon(:,i) = -tau(i)*pid(:,i)/norm(pid(:,i));
    L_tendon(:,i) = -tau(i)*hat(R(:,:,n)*r(:,i))*pid(:,i)/norm(pid(:,i));
end

F_sum = sum(F_tendon,2);
L_sum = sum(L_tendon,2);

F_error = norm(F_sum-n_rod(:,n));
L_error = norm(L_sum-m(:,n));


%Calculate arclength to check solution feasibility
arc = arclength(p(1,:),p(2,:),p(3,:));

%Plot solution
plot3(p(1,:),p(2,:),p(3,:));
xlabel('x');
ylabel('y');
zlabel('z');
grid on
axis([-L,L,-L,L,-L,L]);
title(['Arclength is ',num2str(arc)])