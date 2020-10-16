%Cosserat Model Script based on 2017 paper
%Model multiple backbones with intermediate discs with optimisation
clear variables
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
global r

%Tendon parameters
global r_t
global n_t
global tau

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
rad_s = 0.01; %Radial location of secondary rods from central backbone (m)(approx 10mm)
n = 4; %Number of secondary backbones
r = [0 -rad_s 0 rad_s; rad_s 0 -rad_s 0; 0 0 0 0]; %Radial coordinate profile of secondary backbones through disc (local frame)

%Disc Parameters
nd = 2; %Number of discs (not including base,including end effector)
d = linspace(L/nd,L,nd);  %Disc locations on central backbone

%Tendon Parameters
%Tendon Parameters
n_t = 4; %Number of tendons
rad_t = 0.02; %Radial location of tendons (m)(approx 20mm)

%x,y locations of tendons in rod cross section
r_t = [0 -rad_t 0 rad_t; rad_t 0 -rad_t 0; 0 0 0 0]; %Position vectors of tendons in body frame
%Reference Parameters
%Linear rate of change of frame in reference state
v_ref = [0;0;1];

%% /////////// Model Variables ////////////
%Input force/moments at disc and end effector 
F_end = [0;0;0];
M_end = [0;0;0];

%Input tension
tau = [10 0 0 0]; %Tension for each tendon (N)

%% /////// Initialise Model Variables //////////
%Initial n values are [0;0;0] for all rods at all discs
nm_base = zeros(30,1);

%Initial m values are [0;0;0] for all rods at all discs
nm_disc = zeros(30,1);

%Initial disc intersection based on straight position
s_disc = ones(4,1)*d(1);

%Create initial guess vector (56 elements)
guess = [nm_base;nm_disc;s_disc];
%guess = [2.56960681229459e-06;0.486379053704807;-4.46272496602400;-0.0234049369324981;1.80471070896491e-07;-5.46898103610553e-07;-1.02394018942552e-05;-0.686484351587372;-7.81395347776497;-1.69252647520173e-06;-0.144971043260573;-1.65031655135529;-3.14976671505154e-06;0.489969267764452;5.57726725225122;-2.07797282876364e-06;-0.144974053899226;-1.65027229063539;0.00325929074879909;1.55416405772534e-07;1.10352411584471e-06;-0.0112906482362051;1.03771431624075e-07;-6.23148496823441e-08;-0.0233610985211333;-1.00147708758699e-07;-4.44417317416802e-07;-0.0112901899960900;1.02754422737977e-07;2.01041682031975e-07;-1.89137689148635e-06;-0.187589586400671;-2.02275347985265;0.00995772143090132;-8.78445418140338e-08;4.75426284445245e-07;-7.21130783280786e-06;-0.250662636487879;-2.70254506961868;-5.99675859710489e-07;0.0341351220449601;0.368195380235578;1.06485846808159e-05;0.369976984883335;3.98896217870640;-9.46221309922023e-07;0.0341401133382594;0.368140991313971;0.00837278632542334;-6.12613966753591e-07;-1.00189368432324e-06;0.0142897903934617;-2.13699527481103e-07;8.70868581179099e-08;0.0200647377087361;1.14604273926786e-07;4.30891216597192e-07;0.0142891008843693;-2.59484191851764e-07;-1.72461376900905e-07;0.123385565794129;0.125047030053345;0.126743463642817;0.125047004670198];
%% ///////// Solve Optimisation Problem //////////
%Set fsolve options
options = optimoptions(@fsolve,'Algorithm','levenberg-marquardt','Display','iter-detailed','MaxFunctionEvaluations',50000,'MaxIterations',1000);

%Solve optimisation problem with fsolve
[final_guess,fval,exitflag,output] = fsolve(@MultiShootingMethodTendon,guess,options);

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
plotCircle3D(p_disc,disc_normal',rad_s)
plotCircle3D(pb_L',end_normal',rad_s)

%Graph labels
xlabel('x');
ylabel('y');
zlabel('z');
grid on
axis([-L,L,-L,L,-L,L]);
title(['Arclength is ',num2str(arc)])

%Time stats
time = toc;
fprintf('Algorithm took %.2f minutes to run\n',time/60);