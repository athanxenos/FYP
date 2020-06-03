function [pd,Rd,ud,vd] = f_ode(R,u,v,tau)
%Function for one iteration of ode
%Inputs:
%Systems variables - p,R matrix,v,u 
%Tendon tension - tau vector

%Outputs:
%System variable deriatives - pd,Rd,ud,vd

global K_se
global K_bt
global r
global v_ref
global u_ref
global ud_ref
global vd_ref

n_t = size(r,2);

%Reference position curves for each tendon
pid_b = zeros(3,n_t); %velocity of tendon respect to body frame

%Reset intermediate constants
A = 0;
B = 0;
G = 0;
H = 0;
alpha = 0;
beta =0;

%Iterate through each tendon
for i=1:n_t
    %Tendon curve in body frame
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

c = K_bt*ud_ref-hat(u)*K_bt*(u-u_ref)-hat(v)*K_se*(v-v_ref)-beta;
d = K_se*vd_ref-hat(u)*K_se*(v-v_ref)-alpha;

%Governing system matrix
M = [K_se+A, G; B, K_bt+H];

%Matrix inverse
vu_d = M\[d;c];

%Derivative calculation
vd = vu_d(1:3);
ud = vu_d(4:6);
pd = R*v;
Rd = R*hat(u);

end

