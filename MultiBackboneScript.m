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
%v_total = [v_init;v_init;zeros(16,1)];
v_total = repmat(v_init,10,1);

%Initial U values for all rods at all discs
%Ordered as[backbone_base;backbone_disc;secondary_rods_base;secondary_rods_disc]
%Determine initial shape of rod (C or S)

%u_total = [u_init;-u_init;repmat(u_init,4,1);repmat(-u_init,4,1)];
%u_total = [repmat(u_init,10,1)];
u_total = [u_init;-u_init;repmat(u_init(1:2),4,1);repmat(-u_init(1:2),4,1)];

%Create initial guess vector (52 elements)
guess = [v_total;u_total];

%Result from solved problem
%guess =[-1.00755006832693e-08;4.00814081310874e-08;0.999999995244212;-8.25242544838375e-09;4.53863946609913e-08;1.00000000464015;-7.43460273702034e-08;2.14221732274523e-09;1.35473951072219e-08;-1.23311575844085e-07;2.20805077770712e-07;7.75209218905506e-08;-1.49727652290962e-07;2.05250132458315e-07;2.26389659953645e-07;2.50848869465847e-08;1.20502731597623e-08;1.97282609997371e-07;-1.15948598107838e-07;7.22420833022640e-08;-1.37706604398866e-07;-1.08696249236117e-07;-0.999481843983380;0.000587467052507817;0.0863961868030251;0.999715007181363;-0.0101760455443631;-0.0846995875892702;-0.998637899368204;-0.0628871930475231;-0.668121664461309;-0.937144744202072;0.000659805646877749;-0.0932921166312787;-1.00024631544513;0.0600942489038087;-0.480456319431070;-1.06189428766348;0.00207461873736843;1.08010823461537;0.992824452965230;0.142435120391144;0.652900811844146;0.942130810027191;0.0123924483983522;0.0949875787387578;1.00369737565273;-0.00174078985946176;0.494144210033502;1.05577839310766;-0.143732038361578;-1.07841012935776];
%[u_back_base, u_back_disc, u_sec_base, u_sec_disc] = guess_extract(guess)

%Extract v,u sections from vector
v_guess = guess(1:30);
u_guess = guess(31:end);

%Extract backbone values at base and disc
vb0 = v_guess(1:3);
vbd = v_guess(4:6);
ub0 = u_guess(1:3);
ubd = u_guess(4:6);

%Setup secondary rod matrices
us0 = zeros(3,n);
usd = zeros(3,n);

%Extract v,u values at base and disc for secondary rods
vs0 = reshape(v_guess(7:18),[3,n]);
us0(1:2,:) = reshape(u_guess(7:14),[2,n]);
vsd = reshape(v_guess(19:30),[3,n]);
usd(1:2,:) = reshape(u_guess(15:22),[2,n]);

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

%Calculate n,m at disc for central backbone (global)
nb_d = R_disc*K_se*(vb(end,:)'-v_ref);
mb_d = R_disc*K_bt*ub(end,:)';

%Calculate n,m guesses at first disc for central backbone (global)
nb_dguess = R_disc*K_se*(vbd-v_ref);
mb_dguess = R_disc*K_bt*ubd;

%Initialise constraint/error variables
pos_d = zeros(3,n);
ori_d = zeros(3,n);
n_d = zeros(3,n);
m_d = zeros(3,n);
md_sum = zeros(3,1);
nd_guess = zeros(3,n);
md_guess = zeros(3,n);
E1 = zeros(2,n);
E2 = zeros(2,n);

%Integrate secondary rods until they intersect first disc
for i=1:n
    
    %Find disc intersection points for each rod
    [s_disc(i)] = DiscIntersect(r(:,i),Rs0,vs0(:,i),us0(:,i));
    
    %Integrate secondary rods to first disc
    [ps{i},Rs{i},vs{i},us{i},s_s{i}] = RodODE_Eval(r(:,i),Rs0,vs0(:,i),us0(:,i),0,s_disc(i));
    
    %Calculate positional constraint at disc (local)
    pos_d(:,i) = R_disc'*(ps{i}(end,:)-p_disc)'-r(:,i);
    
    %Positonal Error
    E1(:,i) = pos_d(1:2,i);
    
    %Calculate orientation constraint at disc
    ori_d(:,i) = inv_hat(logm(Rs{i}(:,:,end)'*R_disc));
    
    %Orientation Error
    E2(:,i) = ori_d(1:2,i);
    
    %Calculate n,m at disc (global)
    n_d(:,i) = Rs{i}(:,:,end)*K_se*(vs{i}(end,:)'-v_ref);
    m_d(:,i) = Rs{i}(:,:,end)*K_bt*us{i}(end,:)';
    
    %Calculate n,m guesses at first disc (global)
    nd_guess(:,i) = Rs{i}(:,:,end)*K_se*(vsd(:,i)-v_ref);
    md_guess(:,i) = Rs{i}(:,:,end)*K_bt*usd(:,i);

    %Sum moments at disc
    md_sum = md_sum + cross(ps{i}(end,:)',(nd_guess(:,i) - n_d(:,i))) + md_guess(:,i) - m_d(:,i);
end

%Sum forces at disc end
nd_minus = sum(n_d,2) + nb_d;

%Add guesses for forces at disc start
nd_plus = sum(nd_guess,2) + nb_dguess;

%Force Equilibrium Error
E3 = nd_plus - nd_minus - F_disc;

%Moment Equilibrium Error
E4 = md_sum + cross(p_disc',(nb_dguess - nb_d - F_disc)) + mb_dguess - mb_d - M_disc;

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
    [ps_L,Rs_L,vs_L,us_L,s_c] = RodODE_Eval(ps{i}(end,:)',Rs{i}(:,:,end),vsd(:,i),usd(:,i),s_disc(i),d(2)); 
    
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
    E5(:,i) = Rb_L'*(ps{i}(end,:)' - pb_L)-r(:,i);
    
    %Orientation Error at End Effector
    E6(:,i) = inv_hat(logm(Rs{i}(:,:,end)'*Rb_L));
end

%Force Equilibrium Error
E7 = sum(n_L,2) + nb_L - F_end;

%Moment Equilibrium Error
E8 = mL_sum + cross(pb_L,nb_L) + mb_L - cross(pb_L,F_end) - M_end;

%Combine Residual Vector
residual = [E1(:);E2(:);E3;E4;E5(:);E6(:);E7;E8];
error = norm(residual)
%[PosDisc,OriDisc,FDisc,MDisc,PosEnd,OriEnd,FEnd,MEnd] = residual_extract(residual);
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