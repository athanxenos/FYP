function [yd] = f_ode45(s,y)
%Function to evaluate ODE system at one step for backbone with tendons

%Inputs:
%Systems variables described as state vector - y (p,R,v,u)
%Rod length, not used in equations, used for ode45 input - s

%Outputs:
%State vector derivative - yd (pd,Rd,vd,ud)

%Global variables determined by robot structure
global K_se
global K_bt
global r
global v_ref
global tau

%Unpack state vector
R = reshape(y(4:12),3,3);
v = y(13:15);
u = y(16:18);

%Number of tendons
n_t = 4;

%Initialise key variables
pid_b = zeros(3,n_t); %Velocity of tendon curve respect to body frame

%Reset intermediate constants
A = 0;
B = 0;
G = 0;
H = 0;
alpha = 0;
beta =0;

%Iterate through each tendon

for i=1:n_t
    
    %Tendon curve velocity in body frame
    pid_b(:,i) = hat(u)*r(:,i)+v; 
    
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

c = -hat(u)*K_bt*u-hat(v)*K_se*(v-v_ref)-beta;
d = -hat(u)*K_se*(v-v_ref)-alpha;

%Governing system matrix
M = [K_se+A, G; B, K_bt+H];

%Matrix inverse
vu_d = M\[d;c];

%Derivative calculation
pd = R*v;
Rd = R*hat(u);

%Extract values from vector
vd = vu_d(1:3);
ud = vu_d(4:6);

yd = [pd ; reshape(Rd,9,1); vd; ud];
end
