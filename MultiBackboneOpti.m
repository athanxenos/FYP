%Cosserat Model Script based on 2017 paper
%Model multiple backbones with intermediate discs with optimisation
clear variables
clear global
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
global nd
global n_mid
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
L = 0.285; %Arclength of all rods (m)
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
rad_s = 0.015; %Radial location of secondary rods from central backbone (m)(approx 10mm)
n = 4; %Number of secondary backbones
r = [0 -rad_s 0 rad_s; rad_s 0 -rad_s 0; 0 0 0 0]; %Radial coordinate profile of secondary backbones through disc (local frame)

%Disc Parameters
nd = 7; %Number of discs (not including base or end effector)
n_mid = 4;
d = linspace(0,L,nd+2);  %Disc locations on central backbone

%Reference Parameters
%Linear rate of change of frame in reference state
v_ref = [0;0;1];

%% /////////// Model Variables ////////////
%Input force/moments at disc and end effector 
F_end = [0;0;0];
M_end = [0.2;0;0];
F_disc = [0;0;0];
M_disc = [-0.4;0;0];

%% /////// Initialise Model Variables //////////

nm_guess = zeros((nd+1)*30,1);

%Initial disc intersection based on straight position
s_disc = repmat(d(2:end-1),4,1);

%Create initial guess vector (56 elements)
guess = [nm_guess;s_disc(:)];

%% ///////// Solve Optimisation Problem //////////
%Set fsolve options
options = optimoptions(@fsolve,'Display','iter-detailed','MaxFunctionEvaluations',1000000,'MaxIterations',10000);

%Solve optimisation problem with fsolve
[final_guess,fval,exitflag,output] = fsolve(@MultiShootingMethod,guess,options);

%Run one iteration of code
%[residual] = MultiShootingMethod(guess)
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
for i =1:nd
    plotCircle3D(p_disc(i+1,:),disc_normal(i,:),rad_s)
end

plotCircle3D(pb_L',end_normal',rad_s)

%Graph labels
xlabel('x');
ylabel('y');
zlabel('z');
grid on
axis([-L,L,-L,L,0,L]);
%title(['Arclength is ',num2str(arc)])
title(['M disc is ',num2str(M_disc(1)),'Nm, M end is ',num2str(M_end(1)),'Nm'])

%Time stats
time = toc;
fprintf('Algorithm took %.2f minutes to run\n',time/60);