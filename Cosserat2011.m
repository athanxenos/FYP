%Cosserat Model
clear all
close all
clc 

%Rod Parameters
L = 1; %Arclength of rod (m)
ds = 0.01; %Step size
s = 0:ds:L; %arclength parameter (m)
n = length(s);
r = 0.05; %Radius of rod (m)

%Material parameters
Area = pi*r^2; %Area of cross section (m^2)
EY = 210*10^9; %Young's Modulus (Pa)
Pois = 0.3125; %Poisson ratio (v)
GY = 80*10^9;  %Shear Modulus (Pa) G=E/2(1+v)

%Second moments of area of cross section
Ixx = pi*r^4/4;
Iyy = pi*r^4/4;
Izz = Ixx + Iyy; %Polar moment of inertia

%Stiffness Matrices
K_se = diag([GY*Area GY*Area EY*Area]); %Shear/extension stiffness matrix
K_bt = diag([EY*Ixx EY*Iyy EY*Izz]); %Bending/torsion stiffness matrix

%Reference Parameters
%Linear/angular rate of change of frame in reference state
v_ref = [0;0;1];
u_ref = [0;0;0];
ud_ref = [0;0;0];
vd_ref = [0;0;0];

%Initial Conditions
R0 = eye(3); %Initial rod orientation at base
p0 = [0;0;0]; %Inital rod position at base
%Guess initial conditions for v,u
v0=[0;0;1];
u0=[1;0;0];

%Tendon Parameters
n_t = 4; %Number of tendons
r_t = 0.03; %Radial location of tendons (m)

%x,y locations of tendons in rod cross section
x_t = [0 -r_t 0 r_t];
y_t = [r_t 0 -r_t 0];
r = [x_t; y_t; 0 0 0 0]; %position vectors of tendons

%Reference position curves for each tendon
pi = zeros(3,n_t);
pi_d = zeros(3,n_t); %velocity of tendon respect to body frame

%Tension Input
tau = [1 0 0 0]; %Tension for each tendon


%Setup initial iteration
R=zeros(3,3,n);
p=zeros(3,n);
v=zeros(3,n);
u=zeros(3,n);

%Assign initial values
R(:,:,1)=R0;
p(:,1)=p0;
v(:,1)=v0;
u(:,1)=u0;


A = 0;
B = 0;
G = 0;
H = 0;
alpha = 0;
beta =0;

%n = R*K_se*(v-v_ref);
%m = R*K_bt*(u-u_ref);

for j=1:n-1
    %Tendon path curves and variables (function?)
    for i=1:n_t
    pi(:,i) = R(:,:,j)*r(:,i) + p(:,j);
    pi_d(:,i) = hat(u(:,j))*r(:,i)+v(:,j); %Tendon curve in body frame

    A_i = -tau(i)*(hat(pi_d(:,i))^2)/(norm(pi_d(:,i)))^3;
    A = A + A_i;

    B_i = hat(r(:,i))*A_i;
    B = B + B_i;

    G_i = -A_i*hat(r(:,i));
    G = G + G_i;

    H_i = -B_i*hat(r(:,i));
    H = H + H_i;

    alpha_i = A_i*(hat(u(:,j))*pi_d(:,i));
    alpha = alpha + alpha_i;

    beta_i = hat(r(:,i))*alpha_i;
    beta = beta + beta_i;
    end

    c = K_bt*ud_ref-hat(u(:,j))*K_bt*(u(:,j)-u_ref)-hat(v(:,j))*K_se*(v(:,j)-v_ref)-beta;
    d = K_se*vd_ref-hat(u(:,j))*K_se*(v(:,j)-v_ref)-alpha;

    M = [K_se+A, G; B, K_bt+H];

    vu_d = M\[d;c];

    vd = vu_d(1:3);
    ud = vu_d(4:6);
    pd = R(:,:,j)*v(:,j);
    Rd = R(:,:,j)*hat(u(:,j));

    R(:,:,j+1) = R(:,:,j) + Rd*ds;
    p(:,j+1) = p(:,j) + pd*ds;
    v(:,j+1) = v(:,j) + vd*ds;
    u(:,j+1) = u(:,j) + ud*ds;
end
arclength(p(1,:),p(2,:),p(3,:))
plot3(p(1,:),p(2,:),p(3,:));
grid on