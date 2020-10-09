%Cosserat Model Script based on 2017 paper
%Model multiple backbones with intermediate discs for given state
clear variables
close all
clc 

%Define global variables for model
%Tendon parameters
global r_t
global tau

%ODE parameters
global K_se
global K_bt
global v_ref

%DiscIntersect variables
global d
global R_disc
global p_disc

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

%Tendon Parameters
%Tendon Parameters
n_t = 4; %Number of tendons
rad_t = 0.01; %Radial location of tendons (m)(approx 10mm)

%x,y locations of tendons in rod cross section
r_t = [0 -rad_t 0 rad_t; rad_t 0 -rad_t 0; 0 0 0 0]; %Position vectors of tendons in body frame

%Reference Parameters
%Linear rate of change of frame in reference state
v_ref = [0;0;1];

%% /////////// Model Inputs ////////////
%Input force/moments at disc and end effector 
F_end = [0;0;0];
M_end = [0;0;0];

%Input tension
tau = [0 0 0 0]; %Tension for each tendon (N)

%Set initial n,m values for all rods
nm_base = zeros(30,1);
nm_disc = zeros(30,1);

%% /////// Initialise Model Variables //////////
%Extract backbone values at base and disc
nb0 = nm_base(1:3);
mb0 = nm_base(4:6);
nb_Dguess = nm_disc(1:3);
mb_Dguess = nm_disc(4:6);

%Extract v,u values at base and disc for secondary rods
ns0 = reshape(nm_base(7:18),[3,n]);
ms0 = reshape(nm_base(19:30),[3,n]);
ns_Dguess = reshape(nm_disc(7:18),[3,n]);
ms_Dguess = reshape(nm_disc(19:30),[3,n]);

%Define initial conditions for central backbone
pb0 = [0;0;0];
Rb0 = eye(3);

%Initialise variables for secondary backbones
Rs0 = eye(3);

%Initialise disc intersection arclengths
s_disc = ones(1,4)*d(1);

%Initialise cells for secondary backbones
ps = cell(1,n);
us = cell(1,n);
vs = cell(1,n);
Rs = cell(1,n);
s_s = cell(1,n);

%% /////////// Algorithm Starts ///////////////
%/////////// Base to First Disc /////////////
%Integrate central backbone upto first disc
vb0 = K_se^-1*Rb0'*nb0 + v_ref;
ub0 = K_bt^-1*Rb0'*mb0;

[pb,Rb,vb,ub,s] = TendonODE_Eval(pb0,Rb0,vb0,ub0,0,d(1));

%Extract normal to disc as z axis of backbone R matrix
R_disc = Rb(:,:,end);
p_disc = pb(end,:);
disc_normal = R_disc(:,3);

%Calculate n,m at disc for central backbone (global)
nb_D = R_disc*K_se*(vb(end,:)'-v_ref);
mb_D = R_disc*K_bt*ub(end,:)';

%Initialise key variables
pid = zeros(3,n_t);
F_tendon = zeros(3,n_t);
L_tendon = zeros(3,n_t);

%Iterate through each tendon to calculate boundary force/moment
for i=1:n_t
    pid(:,i) = R_disc*(hat(ub(end,:)')*r_t(:,i)+vb(end,:)');
    F_tendon(:,i) = -tau(i)*pid(:,i)/norm(pid(:,i));
    L_tendon(:,i) = -tau(i)*hat(R_disc*r_t(:,i))*pid(:,i)/norm(pid(:,i));
end

%Sum force/moments
F_disc = sum(F_tendon,2);
M_disc = sum(L_tendon,2);

%Initialise constraint/error variables
pos_D = zeros(3,n);
ori_D = zeros(3,n);
n_D = zeros(3,n);
m_D = zeros(3,n);
mD_sum = zeros(3,1);
E1 = zeros(2,n);
E2 = zeros(2,n);
vs0 = zeros(3,n);
us0 = zeros(3,n);


%Integrate secondary rods until they intersect first disc
for i=1:n
    
    %Find disc intersection points for each rod
    vs0(:,i) = K_se^-1*Rs0'*ns0(:,i) + v_ref;
    us0(:,i) = K_bt^-1*Rs0'*ms0(:,i);
    [s_disc(i)] = DiscIntersect(r(:,i),Rs0,vs0(:,i),us0(:,i));
    
    %Integrate secondary rods to first disc
    [ps{i},Rs{i},vs{i},us{i},s_s{i}] = RodODE_Eval(r(:,i),Rs0,vs0(:,i),us0(:,i),0,s_disc(i));
    
    %Calculate positional constraint at disc (local)
    pos_D(:,i) = R_disc'*(ps{i}(end,:)-p_disc)' - r(:,i);
    
    %Positonal Error
    E1(:,i) = pos_D(1:2,i);
    
    %Calculate orientation constraint at disc
    ori_D(:,i) = inv_hat(logm(Rs{i}(:,:,end)'*R_disc));
    
    %Orientation Error
    E2(:,i) = ori_D(1:2,i);
    
    %Calculate n,m at disc (global)
    n_D(:,i) = Rs{i}(:,:,end)*K_se*(vs{i}(end,:)'-v_ref);
    m_D(:,i) = Rs{i}(:,:,end)*K_bt*us{i}(end,:)';
    
    %Sum moments at disc
    mD_sum = mD_sum + cross(ps{i}(end,:)',(-ns_Dguess(:,i) + n_D(:,i))) - ms_Dguess(:,i) + m_D(:,i);
end

%Sum forces exiting disc
nd_minus = sum(n_D,2) + nb_D;

%Sum forces entering disc
nd_plus = sum(ns_Dguess,2) + nb_Dguess;

%Force Equilibrium Error
E3 = -nd_plus + nd_minus - F_disc;

%Moment Equilibrium Error
E4 = mD_sum + cross(p_disc',(-nb_Dguess + nb_D - F_disc)) - mb_Dguess + mb_D - M_disc;

%% //////////// First Disc to End Effector ////////////
%Integrate central backbone from first disc to end effector
vbD = K_se^-1*Rb(:,:,end)'*nb_Dguess + v_ref;
ubD = K_bt^-1*Rb(:,:,end)'*mb_Dguess;

[pb_end,Rb_end,vb_end,ub_end,s_L] = TendonODE_Eval(pb(end,:)',Rb(:,:,end),vbD,ubD,d(1),d(2));

%Concatenate s,p,R,v,u parameters 
s = [s;s_L(2:end)];
pb = [pb;pb_end(2:end,:)];
Rb = cat(3,Rb,Rb_end(:,:,2:end));
vb = [vb;vb_end(2:end,:)];
ub = [ub;ub_end(2:end,:)];

%Find end effector position/orientation
pb_L = pb(end,:)';       
Rb_L =  Rb(:,:,end);
end_normal = Rb_L(:,3);

%Calculate n,m at end effector for central backbone (global)
nb_L = Rb_L*K_se*(vb(end,:)'-v_ref);
mb_L = Rb_L*K_bt*ub(end,:)';

%Initialise constraint/error variables
n_L = zeros(3,n);
m_L = zeros(3,n);
mL_sum = zeros(3,1);
E5 = zeros(3,n);
E6 = zeros(3,n);
vsD = zeros(3,n);
usD = zeros(3,n);

%Integrate secondary rods from first disc to end effector
for i=1:n
    
    %Integrate secondary rods from first disc to end effector 
    vsD(:,i) = K_se^-1*Rs{i}(:,:,end)'*ns_Dguess(:,i) + v_ref;
    usD(:,i) = K_bt^-1*Rs{i}(:,:,end)'*ms_Dguess(:,i);
    
    [ps_L,Rs_L,vs_L,us_L,s_c] = RodODE_Eval(ps{i}(end,:)',Rs{i}(:,:,end),vsD(:,i),usD(:,i),s_disc(i),d(2)); 
    
    %Concatenate s,p,R,v,u parameters 
    s_s{i} = [s_s{i};s_c(2:end)];
    ps{i} = [ps{i}; ps_L(2:end,:)]; 
    Rs{i} = cat(3,Rs{i},Rs_L(:,:,2:end));
    vs{i} = [vs{i}; vs_L(2:end,:)];
    us{i} = [us{i}; us_L(2:end,:)];
    
    %Calculate n,m values at end effector for secondary rod (global)
    n_L(:,i) = Rs{i}(:,:,end)*K_se*(vs{i}(end,:)'-v_ref);
    m_L(:,i) = Rs{i}(:,:,end)*K_bt*us{i}(end,:)';
    
    %Sum moments at end effector
    mL_sum = mL_sum + cross(ps{i}(end,:)',n_L(:,i)) + m_L(:,i);
    
    %Position Error at end effector (local)
    E5(:,i) = Rb_L'*(ps{i}(end,:)' - pb_L) - r(:,i);
    
    %Orientation Error at End Effector
    E6(:,i) = inv_hat(logm(Rs{i}(:,:,end)'*Rb_L));
end

%Force Equilibrium Error
E7 = sum(n_L,2) + nb_L - F_end;

%Moment Equilibrium Error
E8 = mL_sum + cross(pb_L,nb_L) + mb_L - cross(pb_L,F_end) - M_end;

%Combine Residual Vector (52 length)
residual = [E1(:);E2(:);E3;E4;E5(:);E6(:);E7;E8];

%% ///////////// Plot Solution ////////////////
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