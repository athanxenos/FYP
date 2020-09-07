%Cosserat Model Script based on 2017 paper
%Model multiple backbones
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

%Rod Parameters
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
n_s = 4; %Number of secondary backbones
ps0 = [0 -rad_s 0 rad_s; rad_s 0 -rad_s 0; 0 0 0 0]; %Radial coordinate profile of secondary backbones

%Disc Parameters
d = [L/2;L]; %Disc locations on central backbone (not including base)

%Reference Parameters
%Linear/angular rate of change of frame in reference state
v_ref = [0;0;1];

%/////////// Model Variables ////////////
%Guess initial conditions for v,u for all rods
v_guess=[0;0;1]; %Linear rate of change of frame
u_guess=[-1;0;0]; %Angular rate of change of frame
%////////////////////////////////////////

%Define initial conditions for central backbone
pb0 = [0;0;0];
Rb0 = eye(3);
vb0 = v_guess;
ub0 = u_guess;

%Initialise variables for secondary backbones
Rs0 = eye(3);
vs0 = zeros(3,n_s);
us0 = zeros(3,n_s);

%Initialise disc intersectin arclengths
s_disc = ones(1,4)*d(1);

%Initialise cells for secondary backbones
ps = cell(1,n_s);
us = cell(1,n_s);
vs = cell(1,n_s);
Rs = cell(1,n_s);
s_s = cell(1,n_s);

%/////////// Algorithm Starts////////////////

%Integrate central backbone upto first disc
[pb,Rb,vb,ub,s] = RodODE_Eval(pb0,Rb0,vb0,ub0,0,d(1));

%Extract normal to disc as z axis of backbone R matrix
R_disc = Rb(:,:,end);
p_disc = pb(end,:);
disc_normal = R_disc(:,3);


%Integrate secondary rods until they intersect first disc
for i=1:n_s
    
    %Guess initial conditions at base for each secondary rod
    vs0(:,i) = v_guess;
    us0(:,i) = u_guess;
    
    %Find disc intersection points for each rod
    [s_disc(i)] = DiscIntersect(ps0(:,i),Rs0,vs0(:,i),us0(:,i));
    
    %Integrate secondary rods to first disc
    [ps{i}, Rs{i},vs{i},us{i},s_s{i}] = RodODE_Eval(ps0(:,i),Rs0,vs0(:,i),us0(:,i),0,s_disc(i));

end


%Integrate central backbone from first disc to end effector
[pb_L,Rb_L,vb_L,ub_L,s_L] = RodODE_Eval(pb(end,:)',Rb(:,:,end),vb(end,:)',ub(end,:)',d(1),d(2));

%Concatenate s,p,R,v,u parameters 
s = [s;s_L(2:end)];
pb = [pb;pb_L(2:end,:)];
Rb = cat(3,Rb,Rb_L(:,:,2:end));
vb = [vb;vb_L(2:end,:)];
ub = [ub;ub_L(2:end,:)];

%Integrate secondary rods from first disc to end effector
for i=1:n_s
    
    %Integrate secondary rods from first disc to end effector 
    [ps_L,Rs_L,vs_L,us_L,s_c] = RodODE_Eval(ps{i}(end,:)',Rs{i}(:,:,end),vs{i}(end,:)',us{i}(end,:)',s_disc(i),d(2)); 
    
    %Concatenate s,p,R,v,u parameters 
    s_s{i} = [s_s{i};s_c(2:end)];
    ps{i} = [ps{i}; ps_L(2:end,:)]; 
    Rs{i} = cat(3,Rs{i},Rs_L(:,:,2:end));
    vs{i} = [vs{i}; vs_L(2:end,:)];
    us{i} = [us{i}; us_L(2:end,:)];
end

%Calculate arclength to check solution feasibility
arc = arclength(pb(:,1),pb(:,2),pb(:,3));


centre_index = find(s==L/2);    %Find first disc index
d_centre = pb(centre_index,:);  %Find first disc centre
end_effector = pb(end,:);       %Find end effector position
end_normal = Rb(:,3,end);       %Find end effector normal

%Plot solution for central backbone
hold on
plot3(pb(:,1),pb(:,2),pb(:,3),'b');


%Loop and plot secondary backbones
for i = 1:n_s
    plot3(ps{i}(:,1),ps{i}(:,2),ps{i}(:,3),'r');
end

%Plot first disc and end effector
plotCircle3D(d_centre,disc_normal',rad_s)
plotCircle3D(end_effector,end_normal',rad_s)

%Graph labels
xlabel('x');
ylabel('y');
zlabel('z');
grid on
axis([-L,L,-L,L,-L,L]);
title(['Arclength is ',num2str(arc)])