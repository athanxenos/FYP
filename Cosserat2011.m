%Cosserat Model
clear all
close all
clc 

%Rod Parameters
L = 0.5; %Max arclength of rod (m)
ds = 0.01; %Step size
s = 0:ds:L; %arclength parameter (m)
r = 0.05; %Radius of rod (m)

%Material parameters
Area = pi*r^2; %Area of cross section
EY = 210*10^9; %Young's Modulus (Pa)
Pois = 0.3125; %Poisson ratio (v)
GY = 80*10^9;  %Shear Modulus (Pa) G=E/2(1+v)

%Second moments of area of cross section
Ixx = pi*r^4/4;
Iyy = pi*r^4/4;
Izz = Ixx + Iyy; %Polar moment of inertia

K_se = diag([GY*Area GY*Area EY*Area]); %Shear/extension stiffness matrix
K_bt = diag([EY*Ixx EY*Iyy EY*Izz]); %Bending/torsion stiffness matrix

%Reference Parameters
%Linear/angular rate of change of frame in reference state
v_ref = [0;0;1];
u_ref = [0;0;0];
R_ref = eye(3);
p_ref= @(s) [0;0;s]; %Initial rod position curve


R0 = R_ref; %Initial rod orientation at base
p0 = [0;0;0]; %Inital rod position at base

%Tendon Parameters
n_t = 4; %Number of tendons
r_t = 0.03; %Radial location of tendons (m)

%x,y locations of tendons in rod cross section
x_t = [0 -r_t 0 r_t];
y_t = [r_t 0 -r_t 0];
r = [x_t; y_t; 0 0 0 0]; %position vectors of tendons

%Reference position curves for each tendon
pi_ref = zeros(3,n_t);
pi = zeros(3,n_t);
pi_d = zeros(3,n_t); %velocity of tendon respect to body frame
tau = [1 0 1 0]; %Tension for each tendon

%Start initial iteration


%Initial Conditions
R=R0;
p=p0;
%Guess initial conditions for v,u
v=[0;0;1];
u=[0;0;0];
A = 0;
B = 0;
G = 0;
H = 0;
alpha = 0;
beta =0;

n = R*K_se*(v-v_ref);
m = R*K_bt*(u-u_ref);


%Tendon path curves and variables (function?)
for i=1:n_t
pi_ref(:,i) = R_ref*r(:,i)+ p_ref(0);
pi(:,i) = R*r(:,i) + p;
pi_d(:,i) = hat(u)*r(:,i)+v;

A_i = -tau(i)*(hat(pi_d(:,i))^2)/(norm(pi_d(:,i)))^3;
A = A + A_i;

B_i = hat(r(:,i))*A_i;
B = B + B_i;

G_i = -A_i*hat(r(:,i));
G = G + G_i;

H_i = -B_i*hat(r(:,i));
H = H + H_i;

alpha_i = A_i*(hat(u)*pi_d(:,i));
alpha = alpha + alpha_i;

beta_i = hat(r(:,i))*alpha_i;
beta = beta + beta_i;
end

