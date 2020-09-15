%Cosserat Model Script based on 2017 paper
%Model multiple backbones with intermediate discs
clear all
close all
clc 

%Define global variables for model
global K_se
global K_bt
global v_ref

global d
global R_disc
global p_disc

global pb
global ps

global n
global r
global F_end
global M_end
global F_disc
global M_disc

global disc_normal
global end_normal
global pb_L

global iteration
iteration =0;
tic
%% Rod Parameters
L = 0.25; %Arclength of all rods (m)
rad = 0.0005; %Radius of all rods (m)(approx 0.5mm)
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
K_bt = diag([EY*Ixx EY*Iyy GY*Izz]); %Bending/torsion stiffness matrix

%Rod Parameters
rad_s = 0.015; %Radial location of secondary rods from central backbone (m)(approx 10mm)
n = 4; %Number of secondary backbones
r = [0 -rad_s 0 rad_s; rad_s 0 -rad_s 0; 0 0 0 0]; %Radial coordinate profile of secondary backbones through disc

%Disc Parameters
d = [L/2;L]; %Disc locations on central backbone

%Reference Parameters
%Linear/angular rate of change of frame in reference state
v_ref = [0;0;1];

%% /////////// Model Variables ////////////
%Input force/moments at disc and end effector
F_end = [0;0;0];
M_end = [-0.02;0;0];
F_disc = [0;0;0];
M_disc = [0;0;0];

%Set initial v,u values for all rods
v_init = [0;0;1];
u_init = [0;0;0];

%Set initial v,u conditions for all rods
v_total = [v_init;v_init;zeros(16,1)];
u_total = repmat(u_init,10,1);

%Create initial guess vector
init_guess = [v_total;u_total];

[final_guess,fval,exitflag,output] = fsolve(@MultiShootingMethod,init_guess);
%% ///////////////////////////////////////////////////
%Calculate arclength to check solution feasibility
arc = arclength(pb(:,1),pb(:,2),pb(:,3));

%Plot solution for central backbone
hold on
plot3(pb(:,1),pb(:,2),pb(:,3),'b');

%Loop and plot secondary backbones
for i = 1:n
    plot3(ps{i}(:,1),ps{i}(:,2),ps{i}(:,3),'r');
end

%Plot first disc and end effector
plotCircle3D(p_disc,disc_normal',rad_s)
plotCircle3D(pb_L',end_normal',rad_s)

%Graph labels
xlabel('x');
ylabel('y');
zlabel('z');
grid on
axis([-L,L,-L,L,-L,L]);
title(['Arclength is ',num2str(arc)])
time = toc