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
r = [0 -rad_s 0 rad_s; rad_s 0 -rad_s 0; 0 0 0 0]; %Radial coordinate profile of secondary backbones through disc (local frame)

%Disc Parameters
d = [L/2;L]; %Disc locations on central backbone

%Reference Parameters
%Linear/angular rate of change of frame in reference state
v_ref = [0;0;1];

%% /////////// Model Variables ////////////
%Input force/moments at disc and end effector 
F_end = [0;0;0];
M_end = [0.05;0;0];
F_disc = [0;0;0];
M_disc = [0;0;0];

%Set initial v,u values for all rods
v_init = [0;0;1];
u_init = [-1;0;0];

%Initial V values are [0;0;1] for all rods at all discs
%Secondary rods z element is omitted to keep problem square (52x52)
v_total = [v_init;v_init;zeros(16,1)];

%Initial U values for all rods at all discs
%Ordered as[backbone_base;backbone_disc;secondary_rods_base;secondary_rods_disc]
%Determine initial shape of rod (C or S)
u_total = [u_init;-u_init;repmat(u_init,4,1);repmat(-u_init,4,1)];
%u_total = [repmat(u_init,10,1)];

%Create initial guess vector (52 elements)
init_guess = [v_total;u_total];

%Result from solved problem
%init_guess =[-1.00755006832693e-08;4.00814081310874e-08;0.999999995244212;-8.25242544838375e-09;4.53863946609913e-08;1.00000000464015;-7.43460273702034e-08;2.14221732274523e-09;1.35473951072219e-08;-1.23311575844085e-07;2.20805077770712e-07;7.75209218905506e-08;-1.49727652290962e-07;2.05250132458315e-07;2.26389659953645e-07;2.50848869465847e-08;1.20502731597623e-08;1.97282609997371e-07;-1.15948598107838e-07;7.22420833022640e-08;-1.37706604398866e-07;-1.08696249236117e-07;-0.999481843983380;0.000587467052507817;0.0863961868030251;0.999715007181363;-0.0101760455443631;-0.0846995875892702;-0.998637899368204;-0.0628871930475231;-0.668121664461309;-0.937144744202072;0.000659805646877749;-0.0932921166312787;-1.00024631544513;0.0600942489038087;-0.480456319431070;-1.06189428766348;0.00207461873736843;1.08010823461537;0.992824452965230;0.142435120391144;0.652900811844146;0.942130810027191;0.0123924483983522;0.0949875787387578;1.00369737565273;-0.00174078985946176;0.494144210033502;1.05577839310766;-0.143732038361578;-1.07841012935776];
%[u_back_base, u_back_disc, u_sec_base, u_sec_disc] = guess_extract(init_guess)

options=optimoptions('fsolve','MaxFunctionEvaluations',50000,'MaxIterations',1200);

[final_guess,fval,exitflag,output] = fsolve(@MultiShootingMethod,init_guess,options);
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