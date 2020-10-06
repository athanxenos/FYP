function residual = RodShootingMethod(guess)
%Function that evaluates ODE's for input guess and returns residual when
%compared to known boundary conditions

%Inputs:
% guess - 6x1 vector with initial guess for v,u
% u is initialised as [0;0;0] and v is [0;0;1]

%Outputs:
% residual - returns error in boundary conditions for given guess
% error is 6x1 vector [F_error; L_error]

global K_se
global K_bt
global r
global v_ref
global tau
global p

n_t = 4;
L = 0.25;

%Extarct guess values
v0 = guess(1:3);
u0 = guess(4:6);

%Initial Conditions
R0 = eye(3); %Initial rod orientation at base
p0 = [0;0;0]; %Inital rod position at base

%Initialise key variables
pid = zeros(3,n_t);
F_tendon = zeros(3,n_t);
L_tendon = zeros(3,n_t);

%Setup initial values for ODE
y0 = [p0 ; reshape(R0,9,1); v0; u0];

%Solve ODE from 0 to L
[s,y] = ode45(@f_ode45, [0 L], y0);

%Extract solution curve values
n = length(s);
px = y(:,1);
py = y(:,2);
pz = y(:,3);

%Define curve as global
p= [px py pz];

%Extract state variables at endpoint s=L
R = reshape(y(:,4:12)',3,3,n);
RL = R(:,:,n);
v = [y(:,13) y(:,14) y(:,15)];
vL = v(n,:)';
u = [y(:,16) y(:,17) y(:,18)];
uL = u(n,:)';

%Iterate through each tendon to calculate boundary force/moment
for i=1:n_t
    pid(:,i) = RL*(hat(uL)*r(:,i)+vL);
    F_tendon(:,i) = -tau(i)*pid(:,i)/norm(pid(:,i));
    L_tendon(:,i) = -tau(i)*hat(RL*r(:,i))*pid(:,i)/norm(pid(:,i));
end

%Sum force/moments
F_sum = sum(F_tendon,2);
L_sum = sum(L_tendon,2);

%Evaluate n,m at endpoin
nL = RL*K_se*(vL-v_ref);
mL = RL*K_bt*uL;

%Calculate error between applied force/moment and internal force/moment
F_error = F_sum-nL;
L_error = L_sum-mL;

%Return 6x1 residual
residual = [F_error; L_error];
end

