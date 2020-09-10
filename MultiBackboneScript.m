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
r = [0 -rad_s 0 rad_s; rad_s 0 -rad_s 0; 0 0 0 0]; %Radial coordinate profile of secondary backbones

%Disc Parameters
d = [L/2;L]; %Disc locations on central backbone (not including base)

%Reference Parameters
%Linear/angular rate of change of frame in reference state
v_ref = [0;0;1];

%/////////// Model Variables ////////////
%Guess initial conditions for v,u for all rods
v0_guess = [0;0;1]; %Linear rate of change of frame
u0_guess = [-1;0;0]; %Angular rate of change of frame

vd_guess = [0;0;1];
ud_guess = [-1;0;0];
for i =1:n_s
    nd_guess(:,i) = K_se*(vd_guess-v_ref);
    md_guess(:,i) = K_bt*ud_guess;
end
nb_dguess = K_se*(vd_guess-v_ref);
mb_dguess = K_bt*ud_guess;

F_ext = [0;0;0];
M_ext = [0;0;0];
F_disc = [0;0;0];
M_disc = [0;0;0];
%////////////////////////////////////////

%Define initial conditions for central backbone
pb0 = [0;0;0];
Rb0 = eye(3);
vb0 = v0_guess;
ub0 = u0_guess;

%Initialise variables for secondary backbones
Rs0 = eye(3);
vs0 = zeros(3,n_s);
us0 = zeros(3,n_s);

%Initialise disc intersection arclengths
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

nb_d = R_disc*K_se*(vb(end,:)'-v_ref);
mb_d = R_disc*K_bt*ub(end,:)';

pos_d = zeros(3,n_s);
ori_d = zeros(3,n_s);

n_d = zeros(3,n_s);
m_d = zeros(3,n_s);
md_sum = zeros(3,1);
E1 = zeros(2,n_s);
E2 = zeros(2,n_s);

%Integrate secondary rods until they intersect first disc
for i=1:n_s
    
    %Guess initial conditions at base for each secondary rod
    vs0(:,i) = v0_guess;
    us0(:,i) = u0_guess;
    
    %Find disc intersection points for each rod
    [s_disc(i)] = DiscIntersect(r(:,i),Rs0,vs0(:,i),us0(:,i));
    
    %Integrate secondary rods to first disc
    [ps{i}, Rs{i},vs{i},us{i},s_s{i}] = RodODE_Eval(r(:,i),Rs0,vs0(:,i),us0(:,i),0,s_disc(i));

    pos_d(:,i) = R_disc'*(ps{i}(end,:)-p_disc)'-r(:,i);
    E1(:,i) = pos_d(1:2,i);
    
    ori_d(:,i) = inv_hat(logm(Rs{i}(:,:,end)'*R_disc));
    E2(:,i) = ori_d(1:2,i);
    
    n_d(:,i) = Rs{i}(:,:,end)*K_se*(vs{i}(end,:)'-v_ref);
    m_d(:,i) = Rs{i}(:,:,end)*K_bt*us{i}(end,:)';
    
    md_sum = md_sum + cross(ps{i}(end,:)',(nd_guess(:,i)-n_d(:,i)))+md_guess(:,i)-m_d(:,i);
end
nd_minus = sum(n_d,2)-nb_d;
nd_plus = sum(nd_guess,2)+nb_dguess;
E3 = nd_plus-nd_minus-F_disc;
E4 = md_sum+cross(p_disc',(nb_dguess-nb_d-F_disc))+mb_dguess-mb_d-M_disc;

%Integrate central backbone from first disc to end effector
[pb_L,Rb_L,vb_L,ub_L,s_L] = RodODE_Eval(pb(end,:)',Rb(:,:,end),vb(end,:)',ub(end,:)',d(1),d(2));

%Concatenate s,p,R,v,u parameters 
s = [s;s_L(2:end)];
pb = [pb;pb_L(2:end,:)];
Rb = cat(3,Rb,Rb_L(:,:,2:end));
vb = [vb;vb_L(2:end,:)];
ub = [ub;ub_L(2:end,:)];

d_index = find(s==L/2);    %Find first disc index
d_centre = pb(d_index,:);  %Find first disc centre
pb_L = pb(end,:)';       %Find end effector position
Rb_L =  Rb(:,:,end);
end_normal = Rb_L(:,3);       %Find end effector normal

nb_L = Rb_L*K_se*(vb(end,:)'-v_ref);
mb_L = Rb_L*K_bt*ub(end,:)';

E5 = zeros(3,n_s); %Positional constraint at end effector
E6 = zeros(3,n_s);
n_L = zeros(3,n_s);
m_L = zeros(3,n_s);
mL_sum = zeros(3,1);

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
    
    n_L(:,i) = Rs{i}(:,:,end)*K_se*(vs{i}(end,:)'-v_ref);
    m_L(:,i) = Rs{i}(:,:,end)*K_bt*us{i}(end,:)';
    mL_sum = mL_sum + cross(ps{i}(end,:)',n_L(:,i))+m_L(:,i);
    
    E5(:,i) = Rb_L'*(ps{i}(end,:)'-pb_L)-r(:,i);
    E6(:,i) = inv_hat(logm(Rs{i}(:,:,end)'*Rb_L));
end

E7 = sum(n_L,2)+nb_L-F_ext;
E8 = mL_sum + cross(pb_L,nb_L)+mb_L-cross(pb_L,F_ext)-M_ext;


%Calculate arclength to check solution feasibility
arc = arclength(pb(:,1),pb(:,2),pb(:,3));

%Plot solution for central backbone
hold on
plot3(pb(:,1),pb(:,2),pb(:,3),'b');


%Loop and plot secondary backbones
for i = 1:n_s
    plot3(ps{i}(:,1),ps{i}(:,2),ps{i}(:,3),'r');
end

%Plot first disc and end effector
plotCircle3D(d_centre,disc_normal',rad_s)
plotCircle3D(pb_L',end_normal',rad_s)

%Graph labels
xlabel('x');
ylabel('y');
zlabel('z');
grid on
axis([-L,L,-L,L,-L,L]);
title(['Arclength is ',num2str(arc)])