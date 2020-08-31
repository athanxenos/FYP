%Cosserat Model Script based on 2017 paper
%Model multiple backbones
clear all
close all
clc 

%Define global variables for model
global K_se
global K_bt
global v_ref

%Rod Parameters
L = 0.25; %Arclength of rods (m)
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
K_bt = diag([EY*Ixx EY*Iyy EY*Izz]); %Bending/torsion stiffness matrix

%Rod Parameters
rad_s = 0.015; %Radial location of secondary rods (m)(approx 10mm)
n_s = 4; %Number of secondary backbones
ps0 = [0 -rad_s 0 rad_s; rad_s 0 -rad_s 0; 0 0 0 0]; %Radial coordinate profile of secondary backbones

%Disc Parameters
n_d = 2;
d_L = [L/2;L];

%Reference Parameters
%Linear/angular rate of change of frame in reference state
v_ref = [0;0;1];

%/////////// Model Variables ////////////
%Guess initial conditions for v,u of central backbone
v_guess=[0;0;1]; %Linear rate of change of frame
u_guess=[-1;0;0]; %Angular rate of change of frame
%////////////////////////////////////////

%Define initial conditions for central backbone
pb0 = [0;0;0];
Rb0 = eye(3);
vb0 = v_guess;
ub0 = u_guess;

%Integrate central backbone upto first disc
[pb,Rb,vb,ub] = RodIntegrate(pb0,Rb0,vb0,ub0,d_L(1));


Rs0 = zeros(3,3,n_s);
vs0 = zeros(3,n_s);
us0 = zeros(3,n_s);

%Define initial conditions for secondary backbones
for i=1:n_s
    
    Rs0(:,:,i) = eye(3);
    vs0(:,i) = v_guess;
    us0(:,i) = u_guess;

    [ps_temp,Rs,vs,us] = RodIntegrate(ps0(:,i),Rs0(:,:,i),vs0(:,i),us0(:,i),d_L(1));
    
    ps(:,:,i) = ps_temp;
 
    pxs(:,i) = ps(:,1,i);
    pys(:,i) = ps(:,2,i);
    pzs(:,i) = ps(:,3,i);  
end

%Return solution curve pb for central backbone
px = pb(:,1);
py = pb(:,2);
pz = pb(:,3);

%Calculate arclength to check solution feasibility
arc = arclength(px,py,pz);

%Plot solution
plot3(px,py,pz,pxs,pys,pzs);
xlabel('x');
ylabel('y');
zlabel('z');
grid on
axis([-L,L,-L,L,-L,L]);
title(['Arclength is ',num2str(arc)])