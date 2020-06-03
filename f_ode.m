function [pd,Rd,ud,vd] = f_ode(p0,R0,u0,v0)
%Function for ode
   
%Rod Parameters
rad = 0.0005; %Radius of rod (m)(approx 0.5mm)

%Material parameters
Area = pi*rad^2; %Area of cross section (m^2)
EY = 210*10^9; %Young's Modulus (Pa)
Pois = 0.3125; %Poisson ratio (v)
GY = 80*10^9;  %Shear Modulus (Pa) G=E/2(1+v)

%Second moments of area of cross section
Ixx = pi*rad^4/4;
Iyy = pi*rad^4/4;
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


%Tendon Parameters
n_t = 4; %Number of tendons
r_t = 0.01; %Radial location of tendons (m)(approx 10mm)

%x,y locations of tendons in rod cross section
x_t = [0 -r_t 0 r_t];
y_t = [r_t 0 -r_t 0];
r = [x_t; y_t; 0 0 0 0]; %position vectors of tendons

%Reference position curves for each tendon
pid = zeros(3,n_t);
pid_b = zeros(3,n_t); %velocity of tendon respect to body frame


%Tension Input
tau = [0.1 0 0 0]; %Tension for each tendon

%Assign initial values
R=R0;
p=p0;
v=v0;
u=u0;


A = 0;
B = 0;
G = 0;
H = 0;
alpha = 0;
beta =0;

for i=1:n_t
    pid(:,i) = R*(hat(u)*r(:,i)+v);
    pid_b(:,i) = hat(u)*r(:,i)+v; %Tendon curve in body frame

    A_i = -tau(i)*(hat(pid_b(:,i))^2)/(norm(pid_b(:,i)))^3;
    A = A + A_i;

    B_i = hat(r(:,i))*A_i;
    B = B + B_i;

    G_i = -A_i*hat(r(:,i));
    G = G + G_i;

    H_i = -B_i*hat(r(:,i));
    H = H + H_i;

    alpha_i = A_i*(hat(u)*pid_b(:,i));
    alpha = alpha + alpha_i;

    beta_i = hat(r(:,i))*alpha_i;
    beta = beta + beta_i;
end

c = K_bt*ud_ref-hat(u)*K_bt*(u-u_ref)-hat(v)*K_se*(v-v_ref)-beta;
d = K_se*vd_ref-hat(u)*K_se*(v-v_ref)-alpha;


M = [K_se+A, G; B, K_bt+H];

vu_d = M\[d;c];

vd = vu_d(1:3);
ud = vu_d(4:6);
pd = R*v;
Rd = R*hat(u);

end

