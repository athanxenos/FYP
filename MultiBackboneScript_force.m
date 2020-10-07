%Cosserat Model Script based on 2017 paper
%Model multiple backbones with intermediate discs for given state
clear variables
close all
clc 

%Define global variables for model
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

%Reference Parameters
%Linear rate of change of frame in reference state
v_ref = [0;0;1];

%% /////////// Model Inputs ////////////
%Input force/moments at disc and end effector 
F_end = [0;0;0];
M_end = [0;0;0];
F_disc = [0;0;0];
M_disc = [0;0;0];

%Set initial v,u values for all rods
n_init = [0;0;0];
m_init = [0;0;0];

%% /////// Initialise Model Variables //////////
%Initial n values are [0;0;0] for all rods at all discs
n_total = repmat(n_init,10,1);

%Initial m values are [0;0;0] for all rods at all discs
m_total = repmat(m_init,10,1);

%Create initial guess vector (60 elements)
guess = [n_total;m_total];

%Extract v,u sections from vector
n_guess = guess(1:30);
m_guess = guess(31:end);

%Extract backbone values at base and disc
nb0 = n_guess(1:3);
mb0 = m_guess(1:3);
nb_Dguess = n_guess(4:6);
mb_Dguess = m_guess(4:6);

%Extract v,u values at base and disc for secondary rods
ns0 = reshape(n_guess(7:18),[3,n]);
ms0 = reshape(m_guess(7:18),[3,n]);
ns_Dguess = reshape(n_guess(19:30),[3,n]);
ms_Dguess = reshape(m_guess(19:30),[3,n]);

%Define initial conditions for central backbone
pb0 = [0;0;0];
Rb0 = eye(3);

%Initialise variables for secondary backbones
Rs0 = eye(3);

%Initialise disc intersection arclengths
s_disc = ones(1,4)*d(1);

%Initialise cells for secondary backbones
ps = cell(1,n);
ns = cell(1,n);
ms = cell(1,n);
Rs = cell(1,n);
s_s = cell(1,n);

%% /////////// Algorithm Starts ///////////////
%/////////// Base to First Disc /////////////
%Integrate central backbone upto first disc
[pb,Rb,nb,mb,s] = RodODE_Eval_force(pb0,Rb0,nb0,mb0,0,d(1));

%Extract normal to disc as z axis of backbone R matrix
R_disc = Rb(:,:,end);
p_disc = pb(end,:);
disc_normal = R_disc(:,3);

%Calculate n,m at disc for central backbone (global)
nb_D = nb(end,:)';
mb_D = mb(end,:)';


%Initialise constraint/error variables
pos_D = zeros(3,n);
ori_D = zeros(3,n);
n_D = zeros(3,n);
m_D = zeros(3,n);
mD_sum = zeros(3,1);
E1 = zeros(2,n);
E2 = zeros(2,n);

%Integrate secondary rods until they intersect first disc
for i=1:n
    
    %Find disc intersection points for each rod
    [s_disc(i)] = DiscIntersect_force(r(:,i),Rs0,ns0(:,i),ms0(:,i));
    
    %Integrate secondary rods to first disc
    [ps{i},Rs{i},ns{i},ms{i},s_s{i}] = RodODE_Eval_force(r(:,i),Rs0,ns0(:,i),ms0(:,i),0,s_disc(i));
    
    %Calculate positional constraint at disc (local)
    pos_D(:,i) = R_disc'*(ps{i}(end,:)-p_disc)' - r(:,i);
    
    %Positonal Error
    E1(:,i) = pos_D(1:2,i);
    
    %Calculate orientation constraint at disc
    ori_D(:,i) = inv_hat(logm(Rs{i}(:,:,end)'*R_disc));
    
    %Orientation Error
    E2(:,i) = ori_D(1:2,i);
    
    %Calculate n,m at disc (negative is included as they enter disc)(global)
    n_D(:,i) = ns{i}(end,:)';
    m_D(:,i) = ms{i}(end,:)';
    
    %Sum moments at disc
    mD_sum = mD_sum + cross(ps{i}(end,:)',(ns_Dguess(:,i) + n_D(:,i))) + ms_Dguess(:,i) + m_D(:,i);
end

%Sum forces exiting disc
nd_minus = sum(n_D,2) + nb_D;

%Sum forces entering disc
nd_plus = sum(ns_Dguess,2) + nb_Dguess;

%Force Equilibrium Error
E3 = nd_plus + nd_minus - F_disc;

%Moment Equilibrium Error
E4 = mD_sum + cross(p_disc',(nb_Dguess + nb_D - F_disc)) + mb_Dguess + mb_D - M_disc;

%% //////////// First Disc to End Effector ////////////
%Integrate central backbone from first disc to end effector
[pb_end,Rb_end,nb_end,mb_end,s_L] = RodODE_Eval_force(pb(end,:)',Rb(:,:,end),nb_Dguess,mb_Dguess,d(1),d(2));

%Concatenate s,p,R,v,u parameters 
s = [s;s_L(2:end)];
pb = [pb;pb_end(2:end,:)];
Rb = cat(3,Rb,Rb_end(:,:,2:end));
nb = [nb;nb_end(2:end,:)];
mb = [mb;mb_end(2:end,:)];

%Find end effector position/orientation
pb_L = pb(end,:)';       
Rb_L =  Rb(:,:,end);
end_normal = Rb_L(:,3);

%Calculate n,m at end effector for central backbone (global)
nb_L = nb(end,:)';
mb_L = mb(end,:)';

%Initialise constraint/error variables
n_L = zeros(3,n);
m_L = zeros(3,n);
mL_sum = zeros(3,1);
E5 = zeros(3,n);
E6 = zeros(3,n);

%Integrate secondary rods from first disc to end effector
for i=1:n
    
    %Integrate secondary rods from first disc to end effector 
    [ps_L,Rs_L,ns_L,ms_L,s_c] = RodODE_Eval_force(ps{i}(end,:)',Rs{i}(:,:,end),ns_Dguess(:,i),ms_Dguess(:,i),s_disc(i),d(2)); 
    
    %Concatenate s,p,R,v,u parameters 
    s_s{i} = [s_s{i};s_c(2:end)];
    ps{i} = [ps{i}; ps_L(2:end,:)]; 
    Rs{i} = cat(3,Rs{i},Rs_L(:,:,2:end));
    ns{i} = [ns{i}; ns_L(2:end,:)];
    ms{i} = [ms{i}; ms_L(2:end,:)];
    
    %Calcualte n,m values at end effector for secondary rod (global)
    n_L(:,i) = ns{i}(end,:)';
    m_L(:,i) = ms{i}(end,:)';
    
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