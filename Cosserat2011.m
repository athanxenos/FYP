%Cosserat Model
clear all
close all
clc 

%Rod Parameters
L = 1; %Arclength of rod (m)
ds = 0.01; %Step size
s = 0:ds:L; %arclength parameter (m)
n = length(s);
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

%Initial Conditions
R0 = eye(3); %Initial rod orientation at base
p0 = [0;0;0]; %Inital rod position at base
%Guess initial conditions for v,u
v0=[0;0;1];
u0=[0.5;0;0];

%Tendon Parameters
n_t = 4; %Number of tendons
r_t = 0.01; %Radial location of tendons (m)(approx 10mm)

%x,y locations of tendons in rod cross section
x_t = [0 -r_t 0 r_t];
y_t = [r_t 0 -r_t 0];
r = [x_t; y_t; 0 0 0 0]; %position vectors of tendons

%Reference position curves for each tendon
pid = zeros(3,n_t,n);
pid_b = zeros(3,n_t); %velocity of tendon respect to body frame
F_tendon = zeros(3,n_t);
L_tendon = zeros(3,n_t);

%Tension Input
tau = [0.1 0 0 0]; %Tension for each tendon
L_i = 1; %Tendon termination point
%k=find(s==L_i);
k=n;

%Setup initial iteration
R=zeros(3,3,n);
p=zeros(3,n);
v=zeros(3,n);
u=zeros(3,n);
n_rod=zeros(3,n);
m=zeros(3,n);


%Assign initial values
R(:,:,1)=R0;
p(:,1)=p0;
v(:,1)=v0;
u(:,1)=u0;




for j=1:n
    %Tendon path curves and variables (function?)
    A = 0;
    B = 0;
    G = 0;
    H = 0;
    alpha = 0;
    beta =0;

    for i=1:n_t
    pid(:,i,j) = R(:,:,j)*(hat(u(:,j))*r(:,i)+v(:,j));
    pid_b(:,i) = hat(u(:,j))*r(:,i)+v(:,j); %Tendon curve in body frame

    A_i = -tau(i)*(hat(pid_b(:,i))^2)/(norm(pid_b(:,i)))^3;
    A = A + A_i;

    B_i = hat(r(:,i))*A_i;
    B = B + B_i;

    G_i = -A_i*hat(r(:,i));
    G = G + G_i;

    H_i = -B_i*hat(r(:,i));
    H = H + H_i;

    alpha_i = A_i*(hat(u(:,j))*pid_b(:,i));
    alpha = alpha + alpha_i;

    beta_i = hat(r(:,i))*alpha_i;
    beta = beta + beta_i;
    end

    c = K_bt*ud_ref-hat(u(:,j))*K_bt*(u(:,j)-u_ref)-hat(v(:,j))*K_se*(v(:,j)-v_ref)-beta;
    d = K_se*vd_ref-hat(u(:,j))*K_se*(v(:,j)-v_ref)-alpha;
    
    n_rod(:,j) = R(:,:,j)*K_se*(v(:,j)-v_ref);
    m(:,j) = R(:,:,j)*K_bt*(u(:,j)-u_ref);

    M = [K_se+A, G; B, K_bt+H];

    vu_d = M\[d;c];

    vd = vu_d(1:3);
    ud = vu_d(4:6);
    pd = R(:,:,j)*v(:,j);
    Rd = R(:,:,j)*hat(u(:,j));
    
    if j<102
     R(:,:,j+1) = R(:,:,j) + Rd*ds;
     p(:,j+1) = p(:,j) + pd*ds;
     v(:,j+1) = v(:,j) + vd*ds;
     u(:,j+1) = u(:,j) + ud*ds;
    end
    
   
end

    n_rod(:,j+1) = R(:,:,j+1)*K_se*(v(:,j+1)-v_ref);
    m(:,j+1) = R(:,:,j+1)*K_bt*(u(:,j+1)-u_ref);

for i=1:n_t
    F_tendon(:,i) = -tau(i)*pid(:,i,k)/norm(pid(:,i,k));
    L_tendon(:,i) = -tau(i)*hat(R(:,:,k)*r(:,i))*pid(:,i,k)/norm(pid(:,i,k));
end
F_sum = sum(F_tendon,2);
L_sum = sum(L_tendon,2);

F_err_abs = F_sum-(n_rod(:,k-1)-n_rod(:,k+1));
L_err_abs = L_sum-(m(:,k-1)-m(:,k+1));
F_err_rel = norm(F_err_abs/(n_rod(:,k-1)-n_rod(:,k+1)));
L_err_rel = norm(L_err_abs/(m(:,k-1)-m(:,k+1)));

arclength(p(1,1:end-1),p(2,1:end-1),p(3,1:end-1))

plot3(p(1,:),p(2,:),p(3,:));
xlabel('x');
ylabel('y');
zlabel('z');
grid on
axis([-1,1,-1,1,-1,1]);