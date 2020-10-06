%Cosserat Model Script based on 2017 paper
%Model multiple backbones with intermediate discs for given state
clear variables
close all
clc 

%Define global variables for model
global K_se
global K_bt
global v_ref
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
F_end = [0.5;0;0];
M_end = [0;0;0];
F_disc = [0;0;0];
M_disc = [0;0;0];

%Set initial v,u values for all rods
v_init = [0;0;1];
u_init = [-1;0;0];

%Result from solved problem (F_end = [0.5;0;0])
guess = [-1.73701147009366e-06;-1.12036780868698e-09;0.999997332408388;1.55186603447041e-06;-5.34898089056952e-10;1.00000066400093;-1.74299920479796e-06;-1.24602366546096e-09;0.999997256898494;-2.68766618288772e-06;-6.94466360557488e-09;0.999979287163858;-1.74104670564860e-06;-1.25433755653399e-09;0.999997224962453;-4.90235821163598e-08;1.05653927693667e-08;1.00002889856661;1.60435778687766e-06;-9.62215724871634e-11;0.999999744415863;9.52005640266319e-07;-1.31848505599307e-09;1.00001232587174;1.59920689939696e-06;-1.06279600933628e-10;0.999999777839813;2.22965167562367e-06;2.14630657486568e-09;0.999987706171750;-0.00145483577888189;-0.0549826666963381;0.00100609314430152;0.00179282608243910;0.0198972280434759;-0.000921116733101241;-0.00143749774358323;-0.0564351843038099;-0.000428433486393093;-8.04346078898493e-05;-0.245818249730698;-0.000431751579574590;-0.00143637854898747;-0.0553135395332850;-0.000428129663876874;-0.00325550319887084;0.204848179776368;-0.000374632586186604;0.00176939649493278;0.0293187318229476;0.000513521984517956;0.00148996578485620;-0.0782351430529439;0.000508556908497941;0.00177084836586477;0.0276745011591377;0.000513167986310103;0.00192399387528175;0.147128809549239;0.000467826026556412];

%% /////// Initialise Model Variables //////////
%Initial V values are [0;0;1] for all rods at all discs
v_total = repmat(v_init,10,1);

%Initial U values are [0;0;0] for all rods at all discs
u_total = repmat(u_init,10,1);

%Create initial guess vector (60 elements)
%guess = [v_total;u_total];

%Extract v,u sections from vector
v_guess = guess(1:30);
u_guess = guess(31:end);

%Extract backbone values at base and disc
vb0 = v_guess(1:3);
ub0 = u_guess(1:3);
vbD = v_guess(4:6);
ubD = u_guess(4:6);

%Extract v,u values at base and disc for secondary rods
vs0 = reshape(v_guess(7:18),[3,n]);
us0 = reshape(u_guess(7:18),[3,n]);
vsD = reshape(v_guess(19:30),[3,n]);
usD = reshape(u_guess(19:30),[3,n]);

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

%Calculate n,m at disc for central backbone (global)
nb_D = R_disc*K_se*(vb(end,:)'-v_ref);
mb_D = R_disc*K_bt*ub(end,:)';

%Calculate n,m guesses at first disc for central backbone(negative as they go into disc) (global)
nb_Dguess = R_disc*K_se*(vbD-v_ref);
mb_Dguess = R_disc*K_bt*ubD;

%Initialise constraint/error variables
pos_D = zeros(3,n);
ori_D = zeros(3,n);
n_D = zeros(3,n);
m_D = zeros(3,n);
mD_sum = zeros(3,1);
nD_guess = zeros(3,n);
mD_guess = zeros(3,n);
E1 = zeros(2,n);
E2 = zeros(2,n);

%Integrate secondary rods until they intersect first disc
for i=1:n
    
    %Find disc intersection points for each rod
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
    
    %Calculate n,m at disc (negative is included as they enter disc)(global)
    n_D(:,i) = Rs{i}(:,:,end)*K_se*(vs{i}(end,:)'-v_ref);
    m_D(:,i) = Rs{i}(:,:,end)*K_bt*us{i}(end,:)';
    
    %Calculate n,m guesses at first disc (global)
    nD_guess(:,i) = Rs{i}(:,:,end)*K_se*(vsD(:,i)-v_ref);
    mD_guess(:,i) = Rs{i}(:,:,end)*K_bt*usD(:,i);

    %Sum moments at disc
    mD_sum = mD_sum + cross(ps{i}(end,:)',(nD_guess(:,i) + n_D(:,i))) + mD_guess(:,i) + m_D(:,i);
end

%Sum forces exiting disc
nd_minus = sum(n_D,2) + nb_D;

%Sum forces entering disc
nd_plus = sum(nD_guess,2) + nb_Dguess;

%Force Equilibrium Error
E3 = nd_plus + nd_minus - F_disc;

%Moment Equilibrium Error
E4 = mD_sum + cross(p_disc',(nb_Dguess + nb_D - F_disc)) + mb_Dguess + mb_D - M_disc;

%% //////////// First Disc to End Effector ////////////
%Integrate central backbone from first disc to end effector
[pb_end,Rb_end,vb_end,ub_end,s_L] = RodODE_Eval(pb(end,:)',Rb(:,:,end),vbD,ubD,d(1),d(2));

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

%Integrate secondary rods from first disc to end effector
for i=1:n
    
    %Integrate secondary rods from first disc to end effector 
    [ps_L,Rs_L,vs_L,us_L,s_c] = RodODE_Eval(ps{i}(end,:)',Rs{i}(:,:,end),vsD(:,i),usD(:,i),s_disc(i),d(2)); 
    
    %Concatenate s,p,R,v,u parameters 
    s_s{i} = [s_s{i};s_c(2:end)];
    ps{i} = [ps{i}; ps_L(2:end,:)]; 
    Rs{i} = cat(3,Rs{i},Rs_L(:,:,2:end));
    vs{i} = [vs{i}; vs_L(2:end,:)];
    us{i} = [us{i}; us_L(2:end,:)];
    
    %Calcualte n,m values at end effector for secondary rod (global)
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