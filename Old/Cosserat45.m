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
global tau

%Rod Parameters
L = 0.25; %Arclength of rod (m)
rad = 0.0005; %Radius of central rod (m)(approx 0.5mm)
Area = pi*rad^2; %Area of cross section (m^2)

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

%/////////// Model Variables ////////////
%Guess initial conditions for v,u
v0=[0;0;1]; %Linear rate of change of frame
u0=[0;1;0]; %Angular rate of change of frame

%Tension Input
tau = [0 1 0 0]; %Tension for each tendon (N)
%////////////////////////////////////////

%Initial Conditions
R0 = eye(3); %Initial rod orientation at base
p0 = [0;0;0]; %Inital rod position at base
n0= R0*K_se*(v0-v_ref); %Initial internal force
m0= R0*K_bt*u0; %Initial internal moment

y0 = [p0 ; reshape(R0,9,1); v0; u0];


%Initialise key variables
pid = zeros(3,n_t);
F_tendon = zeros(3,n_t);
L_tendon = zeros(3,n_t);


%Iterate solution from s = [0,L] using ode45

[s,y] = ode45(@f_ode45, [0 L], y0);
n = length(s);
px = y(:,1);
py = y(:,2);
pz = y(:,3);

p= [px py pz];
R = reshape(y(:,4:12)',3,3,n);
RL = R(:,:,n);
v = [y(:,13) y(:,14) y(:,15)];
vL = v(n,:)';
u = [y(:,16) y(:,17) y(:,18)];
uL = u(n,:)';

%Iterate through each tendon
for i=1:n_t
    pid(:,i) = hat(uL)*r(:,i)+vL;
    F_tendon(:,i) = -tau(i)*pid(:,i)/norm(pid(:,i));
    L_tendon(:,i) = -tau(i)*hat(r(:,i))*pid(:,i)/norm(pid(:,i));
end

F_sum = sum(F_tendon,2);
L_sum = sum(L_tendon,2);

%Evaluate n,m at next step (body frame)
nL = K_se*(vL-v_ref);
mL = K_bt*uL

F_error = F_sum-nL;
L_error = L_sum-mL;

residual = [F_error; L_error];


%Calculate arclength to check solution feasibility
arc = arclength(px,py,pz);

%Plot solution
plot3(px,py,pz);
xlabel('x');
ylabel('y');
zlabel('z');
grid on
axis([-L,L,-L,L,-L,L]);
title(['Arclength is ',num2str(arc)])
