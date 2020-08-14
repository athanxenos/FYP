%Cosserat Model Script based on 2011 paper
%Uses ode45 and  fsolve functions to solve model
clear all
close all
clc 

%Define global variables for model
global K_se
global K_bt
global r
global v_ref
global tau
global p

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
u0=[0;0;0]; %Angular rate of change of frame

%Guess input for fsolve
init_guess = [v0; u0];

%Tension Input
tau = [1 0 0 0]; %Tension for each tendon (N)
%////////////////////////////////////////

%Input initial guess into fsolve to optimise solution
final_guess = fsolve(@RodShootingMethod,init_guess)

%Return solution curve p
px = p(:,1);
py = p(:,2);
pz = p(:,3);

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