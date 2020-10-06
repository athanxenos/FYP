%Cosserat Model Script based on 2011 paper
%Plots final robot pose for different input tensions
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

%Guess initial conditions for v,u
v0=[0;0;1]; %Linear rate of change of frame
u0=[0;0;0]; %Angular rate of change of frame

%Guess input for fsolve
init_guess = [v0; u0];

%Define tau step size
step=0:0.5:10;
n=length(step);
tau_range = zeros(n,4);
tau_range(:,1) = step;

for i=1:n
    %Vary Tension Input
    tau = tau_range(i,:); %Tension for each tendon (N)
   

    %Input initial guess into fsolve to optimise solution
    final_guess = fsolve(@RodShootingMethod,init_guess);
    
    %Plot curve
    hold on
    plot3(p(:,1),p(:,2),p(:,3));
end

%Plot solution
xlabel('x');
ylabel('y');
zlabel('z');
grid on
axis([-L,L,-L,L,-L,L]);
title('Robot poses with varying tension (0-10N)')
