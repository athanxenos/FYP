%Cosserat Model Script based on 2017 paper
%Model secondary backbones 
clear all
close all
clc 

%Define global variables for model
global K_se
global K_bt
global v_ref
global p
global L_L


%Rod Parameters
L = 0.25; %Arclength of rod (m)
rad = 0.0005; %Radius of secondary rod (m)(approx 0.5mm)
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

%Rod Parameters
r_t = 0.015; %Radial location of secondary rod (m)(approx 10mm)


%Reference Parameters
%Linear/angular rate of change of frame in reference state
v_ref = [0;0;1];

%/////////// Model Variables ////////////
%Guess initial conditions for v,u
v0=[0;0;1]; %Linear rate of change of frame
u0=[0;0;0]; %Angular rate of change of frame

L_L = [-0.01;0;0];
%////////////////////////////////////////

%Guess input for fsolve
init_guess = [v0; u0];


%////////////////////////////////////////

%Input initial guess into fsolve to optimise solution
final_guess = fsolve(@RodShootingMethodSec,init_guess);

v_final = final_guess(1:3)
u_final = final_guess(4:6)

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