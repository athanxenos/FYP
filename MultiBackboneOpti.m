%Cosserat Model Script based on 2017 paper
%Model multiple backbones with intermediate discs with optimisation
clear variables
close all
clc 

%Define global variables for model
%ODE parameters
global K_se
global K_bt
global v_ref

%Model Parameters
global d
global n
global r

%Input Variables
global F_end
global M_end
global F_disc
global M_disc

%Plotting Variables
global pb
global ps
global p_disc
global pb_L
global disc_normal
global end_normal

%Start timer
tic
%% //////////// Model Parameters ////////////
%Rod Parameters
L = 0.25; %Arclength of all rods (m)
rad = 0.0005; %Radius of all rods (m)(approx 0.5mm)
Area = pi*rad^2; %Area of cross section (m^2)

%Material Parameters
EY = 210*10^9; %Young's Modulus (Pa)
GY = 80*10^9;  %Shear Modulus (Pa) G=E/2(1+v) (v=0.3125)

%Second moments of area of cross section
Ixx = pi*rad^4/4;
Iyy = pi*rad^4/4;
Izz = Ixx + Iyy; %Polar moment of inertia

%Stiffness Matrices
K_se = diag([GY*Area GY*Area EY*Area]); %Shear/extension stiffness matrix
K_bt = diag([EY*Ixx EY*Iyy GY*Izz]); %Bending/torsion stiffness matrix

%Secondary Rod Parameters
rad_s = 0.015; %Radial location of secondary rods from central backbone (m)(approx 15mm)
n = 4; %Number of secondary backbones
r = [0 -rad_s 0 rad_s; rad_s 0 -rad_s 0; 0 0 0 0]; %Radial coordinate profile of secondary backbones through disc (local frame)

%Disc Parameters
nd = 2; %Number of discs (not including base,including end effector)
d = linspace(L/nd,L,nd);  %Disc locations on central backbone

%Reference Parameters
%Linear rate of change of frame in reference state
v_ref = [0;0;1];

%% /////////// Model Variables ////////////
%Input force/moments at disc and end effector 
F_end = [0.5;0;0];
M_end = [0;0;0];
F_disc = [0;0;0];
M_disc = [0;0;0];

%Set initial v,u values for all rods
v_init = [0;0;1];
u_init = [0;0;0];

%% /////// Initialise Model Variables //////////
%Initial n values are [0;0;0] for all rods at all discs
nm_base = zeros(30,1);

%Initial m values are [0;0;0] for all rods at all discs
nm_disc = zeros(30,1);

%Initial disc intersection based on straight position
s_disc = ones(4,1)*d(1);

%Create initial guess vector (56 elements)
guess = [nm_base;nm_disc;s_disc];

%% ///////// Solve Optimisation Problem //////////
%Set fsolve options
options = optimoptions(@fsolve,'Algorithm','levenberg-marquardt','Display','iter-detailed','MaxFunctionEvaluations',100000,'MaxIterations',1000);

%Solve optimisation problem with fsolve
[final_guess,fval,exitflag,output] = fsolve(@MultiShootingMethod,guess,options);

%% /////////// Plot Solution //////////////
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

%Time stats
time = toc;
fprintf('Algorithm took %.2f minutes to run\n',time/60);