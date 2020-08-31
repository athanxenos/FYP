function residual = RodShootingMethodSec(guess)
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
global v_ref
global p
global L_L

L = 0.25;
F_L = [0;0;0];

%Extract guess values
v0 = guess(1:3);
u0 = guess(4:6);

%Initial Conditions
R0 = eye(3); %Initial rod orientation at base
p0 = [0;0;0]; %Inital rod position at base

%Setup initial values for ODE
y0 = [p0 ; reshape(R0,9,1); v0; u0];

%Solve ODE from 0 to L
[s,y] = ode45(@f_secondary, [0 L], y0);

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

%Evaluate n,m at endpoint
nL = RL*K_se*(vL-v_ref);
mL = RL*K_bt*uL;

%Calculate error between applied force/moment and internal force/moment
F_error = F_L-nL;
L_error = L_L-mL;

%Return 6x1 residual
residual = [F_error; L_error];
end

