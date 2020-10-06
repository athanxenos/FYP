function [residual] = MultiShootingMethod(guess)
%Function that evaluates Multi Backbone ODE's for input guess and returns
%residual of position/orientation/equilibrium conditions

%Inputs:
% guess - 64x1 vector with initial guess for v,u,s_disc (52 length minimum)
% All v/u values are initialised as [0;0;1] and [0;0;0]
% s_disc intersections are initialised for straight rod position

%Outputs:
% residual - returns error vector for given guess (56 length minimum)

%Define global variables for model
%ODE parameters
global K_se
global K_bt
global v_ref

%Model Parameters
global d
global n
global r

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

%Extract v,u sections from vector
v_guess = guess(1:30);
u_guess = guess(31:end);

%Extract backbone values at base and disc
vb0 = v_guess(1:3);
ub0 = u_guess(1:3);
vbD = v_guess(4:6);
ubD = u_guess(4:6);

%Extract v,u values at base and disc for secondary rods
vs0 = reshape(v_guess(7:18),[3,n]);
us0 = reshape(u_guess(7:18),[3,n]);
vsD = reshape(v_guess(19:30),[3,n]);
usD = reshape(u_guess(19:30),[3,n]);

%Extract disc intersection values
s_disc = guess(end-3:end);

%Define initial conditions for central backbone
pb0 = [0;0;0];
Rb0 = eye(3);

%Initialise variables for secondary backbones
Rs0 = eye(3);

%Initialise cells for secondary backbones
ps = cell(1,n);
us = cell(1,n);
vs = cell(1,n);
Rs = cell(1,n);
s_s = cell(1,n);

%Integrate central backbone upto first disc
[pb,Rb,vb,ub,s] = RodODE_Eval(pb0,Rb0,vb0,ub0,0,d(1));

%Extract normal to disc as z axis of backbone R matrix
R_disc = Rb(:,:,end);
p_disc = pb(end,:);
disc_normal = R_disc(:,3);

%Calculate n,m at disc for central backbone (global)
nb_D = R_disc*K_se*(vb(end,:)'-v_ref);
mb_D = R_disc*K_bt*ub(end,:)';

%Calculate n,m guesses at first disc for central backbone(negative as they go into disc) (global)
nb_Dguess = R_disc*K_se*(vbD-v_ref);
mb_Dguess = R_disc*K_bt*ubD;

%Initialise constraint/error variables
pos_D = zeros(3,n);
ori_D = zeros(3,n);
n_D = zeros(3,n);
m_D = zeros(3,n);
mD_sum = zeros(3,1);
nD_guess = zeros(3,n);
mD_guess = zeros(3,n);
E1 = zeros(2,n);
E2 = zeros(2,n);
E_inter = zeros(1,n);

%Integrate secondary rods until they intersect first disc
for i=1:n
  
    %Integrate secondary rods to first disc
    [ps{i},Rs{i},vs{i},us{i},s_s{i}] = RodODE_Eval(r(:,i),Rs0,vs0(:,i),us0(:,i),0,s_disc(i));
    
    %Calculate disc intersection error
    intersect = R_disc.'*(ps{i}(end,:)-p_disc)';
    E_inter(:,i) = intersect(end);
    
    %Calculate positional constraint at disc (local)
    pos_D(:,i) = R_disc'*(ps{i}(end,:)-p_disc)' - r(:,i);
    
    %Positonal Error
    E1(:,i) = pos_D(1:2,i);
    
    %Calculate orientation constraint at disc
    ori_D(:,i) = inv_hat(logm(Rs{i}(:,:,end)'*R_disc));
    
    %Orientation Error
    E2(:,i) = ori_D(1:2,i);
    
    %Calculate n,m at disc (negative is included as they enter disc)(global)
    n_D(:,i) = Rs{i}(:,:,end)*K_se*(vs{i}(end,:)'-v_ref);
    m_D(:,i) = Rs{i}(:,:,end)*K_bt*us{i}(end,:)';
    
    %Calculate n,m guesses at first disc (global)
    nD_guess(:,i) = Rs{i}(:,:,end)*K_se*(vsD(:,i)-v_ref);
    mD_guess(:,i) = Rs{i}(:,:,end)*K_bt*usD(:,i);

    %Sum moments at disc
    mD_sum = mD_sum + cross(ps{i}(end,:)',(nD_guess(:,i) + n_D(:,i))) + mD_guess(:,i) + m_D(:,i);
end

%Sum forces exiting disc
nd_minus = sum(n_D,2) + nb_D;

%Sum forces entering disc
nd_plus = sum(nD_guess,2) + nb_Dguess;

%Force Equilibrium Error
E3 = nd_plus + nd_minus - F_disc;

%Moment Equilibrium Error
E4 = mD_sum + cross(p_disc',(nb_Dguess + nb_D - F_disc)) + mb_Dguess + mb_D - M_disc;

%Integrate central backbone from first disc to end effector
[pb_end,Rb_end,vb_end,ub_end,s_L] = RodODE_Eval(pb(end,:)',Rb(:,:,end),vbD,ubD,d(1),d(2));

%Concatenate s,p,R,v,u parameters 
s = [s;s_L(2:end)];
pb = [pb;pb_end(2:end,:)];
Rb = cat(3,Rb,Rb_end(:,:,2:end));
vb = [vb;vb_end(2:end,:)];
ub = [ub;ub_end(2:end,:)];

%Find end effector position/orientation
pb_L = pb(end,:)';       
Rb_L =  Rb(:,:,end);
end_normal = Rb_L(:,3);

%Calculate n,m at end effector for central backbone (global)
nb_L = Rb_L*K_se*(vb(end,:)'-v_ref);
mb_L = Rb_L*K_bt*ub(end,:)';

%Initialise constraint/error variables
n_L = zeros(3,n);
m_L = zeros(3,n);
mL_sum = zeros(3,1);
E5 = zeros(3,n);
E6 = zeros(3,n);

%Integrate secondary rods from first disc to end effector
for i=1:n
    
    %Integrate secondary rods from first disc to end effector 
    [ps_L,Rs_L,vs_L,us_L,s_c] = RodODE_Eval(ps{i}(end,:)',Rs{i}(:,:,end),vsD(:,i),usD(:,i),s_disc(i),d(2)); 
    
    %Concatenate s,p,R,v,u parameters 
    s_s{i} = [s_s{i};s_c(2:end)];
    ps{i} = [ps{i}; ps_L(2:end,:)]; 
    Rs{i} = cat(3,Rs{i},Rs_L(:,:,2:end));
    vs{i} = [vs{i}; vs_L(2:end,:)];
    us{i} = [us{i}; us_L(2:end,:)];
    
    %Calcualte n,m values at end effector for secondary rod (global)
    n_L(:,i) = Rs{i}(:,:,end)*K_se*(vs{i}(end,:)'-v_ref);
    m_L(:,i) = Rs{i}(:,:,end)*K_bt*us{i}(end,:)';
    
    %Sum moments at end effector
    mL_sum = mL_sum + cross(ps{i}(end,:)',n_L(:,i)) + m_L(:,i);
    
    %Position Error at end effector (local)
    E5(:,i) = Rb_L'*(ps{i}(end,:)' - pb_L) - r(:,i);
    
    %Orientation Error at End Effector
    E6(:,i) = inv_hat(logm(Rs{i}(:,:,end)'*Rb_L));
end

%Force Equilibrium Error
E7 = sum(n_L,2) + nb_L - F_end;

%Moment Equilibrium Error
E8 = mL_sum + cross(pb_L,nb_L) + mb_L - cross(pb_L,F_end) - M_end;

%Combine Residual Vector
residual = [E_inter(:);E1(:);E2(:);E3;E4;E5(:);E6(:);E7;E8];

end

