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
M_end = [0;0;0];
F_disc = [0;0;0];
M_disc = [0;0;0];

%Set initial v,u values for all rods
v_init = [0;0;1];
u_init = [-1;0;0];

%Set initial v,u conditions for all rods
v_total = [v_init;v_init;zeros(16,1)];
u_total = repmat(u_init,10,1);

%Create initial guess vector
init_guess = [v_total;u_total];

%Extract v,u sections from vector
v_guess = init_guess(1:22);
u_guess = init_guess(23:52);

%Extract backbone values
vb0 = v_guess(1:3);
vbd = v_guess(4:6);
ub0 = u_guess(1:3);
ubd = u_guess(4:6);

%Setup secondary rod matrices
vs0 = ones(3,n);
vsd = ones(3,n);

%Extract v,u values at base and disc
vs0(1:2,:) = reshape(v_guess(7:14),[2,n]);
us0 = reshape(u_guess(7:18),[3,n]);
vsd(1:2,:) = reshape(v_guess(15:end),[2,n]);
usd = reshape(u_guess(19:end),[3,n]);

%Setup n,m matrices
nd_guess = zeros(3,n);
md_guess = zeros(3,n);

%Calculate n,m values at first disc for secondary rods
for i =1:n
    nd_guess(:,i) = K_se*(vsd(:,i)-v_ref);
    md_guess(:,i) = K_bt*usd(:,i);
end

%Calculate n,m values at first disc for central backbone
nb_dguess = K_se*(vbd-v_ref);
mb_dguess = K_bt*ubd;

%% ////////////////////////////////////////

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
[pb,Rb,vb,ub,s] = RodODE_Eval(pb0,Rb0,vb0,ub0,0,d(1));

%Extract normal to disc as z axis of backbone R matrix
R_disc = Rb(:,:,end);
p_disc = pb(end,:);
disc_normal = R_disc(:,3);

%Calculate n,m at disc for central backbone
nb_d = R_disc*K_se*(vb(end,:)'-v_ref);
mb_d = R_disc*K_bt*ub(end,:)';

%Initialise constraint/error variables
pos_d = zeros(3,n);
ori_d = zeros(3,n);
n_d = zeros(3,n);
m_d = zeros(3,n);
md_sum = zeros(3,1);
E1 = zeros(2,n);
E2 = zeros(2,n);

%Integrate secondary rods until they intersect first disc
for i=1:n
    
    %Find disc intersection points for each rod
    [s_disc(i)] = DiscIntersect(r(:,i),Rs0,vs0(:,i),us0(:,i));
    
    %Integrate secondary rods to first disc
    [ps{i},Rs{i},vs{i},us{i},s_s{i}] = RodODE_Eval(r(:,i),Rs0,vs0(:,i),us0(:,i),0,s_disc(i));
    
    %Calculate positional constraint at disc
    pos_d(:,i) = R_disc'*(ps{i}(end,:)-p_disc)'-r(:,i);
    
    %Positonal Error
    E1(:,i) = pos_d(1:2,i);
    
    %Calculate orientation constraint at disc
    ori_d(:,i) = inv_hat(logm(Rs{i}(:,:,end)'*R_disc));
    
    %Orientation Error
    E2(:,i) = ori_d(1:2,i);
    
    %Calculate n,m at disc
    n_d(:,i) = Rs{i}(:,:,end)*K_se*(vs{i}(end,:)'-v_ref);
    m_d(:,i) = Rs{i}(:,:,end)*K_bt*us{i}(end,:)';
    
    %Sum moments at disc
    md_sum = md_sum + cross(ps{i}(end,:)',(nd_guess(:,i)-n_d(:,i)))+md_guess(:,i)-m_d(:,i);
end

%Sum forces at disc end
nd_minus = sum(n_d,2)-nb_d;

%Add guesses for forces at disc start
nd_plus = sum(nd_guess,2)+nb_dguess;

%Force Equilibrium Error
E3 = nd_plus-nd_minus-F_disc;

%Moment Equilibrium Error
E4 = md_sum+cross(p_disc',(nb_dguess-nb_d-F_disc))+mb_dguess-mb_d-M_disc;

%% //////////// First Disc to End Effector ////////////
%Integrate central backbone from first disc to end effector
[pb_L,Rb_L,vb_L,ub_L,s_L] = RodODE_Eval(pb(end,:)',Rb(:,:,end),vbd,ubd,d(1),d(2));

%Concatenate s,p,R,v,u parameters 
s = [s;s_L(2:end)];
pb = [pb;pb_L(2:end,:)];
Rb = cat(3,Rb,Rb_L(:,:,2:end));
vb = [vb;vb_L(2:end,:)];
ub = [ub;ub_L(2:end,:)];

%Find end effector position/orientation
pb_L = pb(end,:)';       
Rb_L =  Rb(:,:,end);
end_normal = Rb_L(:,3);

%Calculate n,m at end effector for central backbone
nb_L = Rb_L*K_se*(vb(end,:)'-v_ref);
mb_L = Rb_L*K_bt*ub(end,:)';

%Initialise constraint/error variables
n_L = zeros(3,n);
m_L = zeros(3,n);
mL_sum = zeros(3,1);
E5 = zeros(3,n);
E6 = zeros(3,n);

%Integrate secondary rods from first disc to end effector
for i=1:n
    
    %Integrate secondary rods from first disc to end effector 
    [ps_L,Rs_L,vs_L,us_L,s_c] = RodODE_Eval(ps{i}(end,:)',Rs{i}(:,:,end),vsd(:,i),usd(:,i),s_disc(i),d(2)); 
    
    %Concatenate s,p,R,v,u parameters 
    s_s{i} = [s_s{i};s_c(2:end)];
    ps{i} = [ps{i}; ps_L(2:end,:)]; 
    Rs{i} = cat(3,Rs{i},Rs_L(:,:,2:end));
    vs{i} = [vs{i}; vs_L(2:end,:)];
    us{i} = [us{i}; us_L(2:end,:)];
    
    %Calcualte n,m values at end effector for secondary rod
    n_L(:,i) = Rs{i}(:,:,end)*K_se*(vs{i}(end,:)'-v_ref);
    m_L(:,i) = Rs{i}(:,:,end)*K_bt*us{i}(end,:)';
    
    %Sum moments at end effector
    mL_sum = mL_sum + cross(ps{i}(end,:)',n_L(:,i))+m_L(:,i);
    
    %Position Error at end effector
    E5(:,i) = Rb_L'*(ps{i}(end,:)'-pb_L)-r(:,i);
    
    %Orientation Error at End Effector
    E6(:,i) = inv_hat(logm(Rs{i}(:,:,end)'*Rb_L));
end

%Force Equilibrium Error
E7 = sum(n_L,2)+nb_L-F_end;

%Moment Equilibrium Error
E8 = mL_sum + cross(pb_L,nb_L)+mb_L-cross(pb_L,F_end)-M_end;

%Combine Residual Vector
residual = [E1(:);E2(:);E3;E4;E5(:);E6(:);E7;E8];

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