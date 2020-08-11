function residual = RodShootingMethod(guess)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
global K_se
global K_bt
global r
global v_ref
global tau
global p
global v0
n_t = 4;
L = 0.25;

%v0 = guess(1:3);
u0 = guess;

%Initial Conditions
R0 = eye(3); %Initial rod orientation at base
p0 = [0;0;0]; %Inital rod position at base

%Initialise key variables
pid = zeros(3,n_t);
F_tendon = zeros(3,n_t);
L_tendon = zeros(3,n_t);


y0 = [p0 ; reshape(R0,9,1); v0; u0];

[s,y] = ode45(@f_ode45, [0 L], y0);

n = length(s);
px = y(:,1);
py = y(:,2);
pz = y(:,3);

p= [px py pz];
R = reshape(y(:,4:12)',3,3,n);
RL = R(:,:,n);
v = [y(:,13) y(:,14) y(:,15)];
vL = v(n,:)';
u = [y(:,16) y(:,17) y(:,18)];
uL = u(n,:)';

%Iterate through each tendon
for i=1:n_t
    pid(:,i) = RL*(hat(uL)*r(:,i)+vL);
    F_tendon(:,i) = -tau(i)*pid(:,i)/norm(pid(:,i));
    L_tendon(:,i) = -tau(i)*hat(RL*r(:,i))*pid(:,i)/norm(pid(:,i));
end

F_sum = sum(F_tendon,2);
L_sum = sum(L_tendon,2);

%Evaluate n,m at next step
nL = RL*K_se*(vL-v_ref);
mL = RL*K_bt*uL;

F_error = norm(F_sum-nL);
L_error = norm(L_sum-mL);

residual = [L_error]
end

