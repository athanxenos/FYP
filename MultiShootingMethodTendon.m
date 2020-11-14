function [residual] = MultiShootingMethodTendon(guess)
%MultiShootingMethod
%Written by Athan Xenos

%Function that evaluates Multi Backbone ODE's for input pose guess and returns
%residual of position/orientation/equilibrium conditions
%Uses coupled tendon method

%Inputs:
% guess - 64x1 vector with initial guess for v,u,s_disc 
% All v/u values are initialised as [0;0;1] and [0;0;0]
% s_disc intersections are initialised for straight rod position

%Outputs:
% residual - returns error vector for given guess 

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

%Tendon parameters
global r_t
global n_t
global tau

%Plotting Variables
global pb
global ps
global p_disc
global pb_L
global disc_normal
global end_normal

%Extract v,u sections from vector
nm_base = guess(1:30);
nm_disc = guess(31:60);

%Extract backbone values at base and disc
nb0 = nm_base(1:3);
mb0 = nm_base(4:6);
nb_Dguess = nm_disc(1:3);
mb_Dguess = nm_disc(4:6);

%Extract v,u values at base and disc for secondary rods
ns0 = reshape(nm_base(7:18),[3,n]);
ms0 = reshape(nm_base(19:30),[3,n]);
ns_Dguess = reshape(nm_disc(7:18),[3,n]);
ms_Dguess = reshape(nm_disc(19:30),[3,n]);

%Extract disc intersection values
s_disc = guess(end-3:end);

%Define initial conditions for central backbone
pb0 = [0;0;0];
Rb0 = eye(3);

%Initialise variables for secondary backbones
Rs0 = eye(3);

%Initialise cells for secondary backbones
ps = cell(1,n);
ms = cell(1,n);
ns = cell(1,n);
Rs = cell(1,n);
s_s = cell(1,n);

%% /////////// Base to First Disc /////////////
%Integrate central backbone upto first disc
vb0 = K_se^-1*Rb0'*nb0 + v_ref;
ub0 = K_bt^-1*Rb0'*mb0;

[pb,Rb,vb,ub,s] = TendonODE_Eval(pb0,Rb0,vb0,ub0,0,d(1));

%Extract normal to disc as z axis of backbone R matrix
R_disc = Rb(:,:,end);
p_disc = pb(end,:);
disc_normal = R_disc(:,3);

%Calculate n,m at disc for central backbone (global)
nb_D = R_disc*K_se*(vb(end,:)'-v_ref);
mb_D = R_disc*K_bt*ub(end,:)';

%Initialise key variables
pid = zeros(3,n_t);
F_tendon = zeros(3,n_t);
L_tendon = zeros(3,n_t);

%Iterate through each tendon to calculate boundary force/moment
for i=1:n_t
    pid(:,i) = R_disc*(hat(ub(end,:)')*r_t(:,i)+vb(end,:)');
    F_tendon(:,i) = -tau(i)*pid(:,i)/norm(pid(:,i));
    L_tendon(:,i) = -tau(i)*hat(R_disc*r_t(:,i))*pid(:,i)/norm(pid(:,i));
end

%Sum force/moments
F_disc = sum(F_tendon,2);
M_disc = sum(L_tendon,2);

%Initialise constraint/error variables
pos_D = zeros(3,n);
ori_D = zeros(3,n);
n_D = zeros(3,n);
m_D = zeros(3,n);
mD_sum = zeros(3,1);
E1 = zeros(2,n);
E2 = zeros(2,n);
E_inter = zeros(1,n);


%Integrate secondary rods until they intersect first disc
for i=1:n
  
    %Integrate secondary rods to first disc
    [ps{i},Rs{i},ns{i},ms{i},s_s{i}] = CosseratODE_Eval(r(:,i),Rs0,ns0(:,i),ms0(:,i),0,s_disc(i));
    
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
    n_D(:,i) = ns{i}(end,:)';
    m_D(:,i) = ms{i}(end,:)';

    %Sum moments at disc
    mD_sum = mD_sum + cross(ps{i}(end,:)',(-ns_Dguess(:,i) + n_D(:,i))) - ms_Dguess(:,i) + m_D(:,i);
end

%Sum forces exiting disc
nd_minus = sum(n_D,2) + nb_D;

%Sum forces entering disc
nd_plus = sum(ns_Dguess,2) + nb_Dguess;

%Force Equilibrium Error
E3 = -nd_plus + nd_minus - F_disc;

%Moment Equilibrium Error
E4 = mD_sum + cross(p_disc',(-nb_Dguess + nb_D - F_disc)) - mb_Dguess + mb_D - M_disc;

%% //////////// First Disc to End Effector ////////////
%Integrate central backbone from first disc to end effector
[pb_end,Rb_end,nb,mb,s_L] = CosseratODE_Eval(pb(end,:)',Rb(:,:,end),nb_Dguess,mb_Dguess,d(1),d(2));

%Concatenate s,p,R,v,u parameters 
s = [s;s_L(2:end)];
pb = [pb;pb_end(2:end,:)];
Rb = cat(3,Rb,Rb_end(:,:,2:end));

%Find end effector position/orientation
pb_L = pb(end,:)';       
Rb_L =  Rb(:,:,end);
end_normal = Rb_L(:,3);

%Calculate n,m at end effector for central backbone (global)
%nb_L = Rb_L*K_se*(vb(end,:)'-v_ref);
%mb_L = Rb_L*K_bt*ub(end,:)';
nb_L = nb(end,:)';
mb_L = mb(end,:)';

%Initialise constraint/error variables
n_L = zeros(3,n);
m_L = zeros(3,n);
mL_sum = zeros(3,1);
E5 = zeros(3,n);
E6 = zeros(3,n);

%Integrate secondary rods from first disc to end effector
for i=1:n
    
    %Integrate secondary rods from first disc to end effector 
    [ps_L,Rs_L,ns_L,ms_L,s_c] = CosseratODE_Eval(ps{i}(end,:)',Rs{i}(:,:,end),ns_Dguess(:,i),ms_Dguess(:,i),s_disc(i),d(2)); 
    
    %Concatenate s,p,R,v,u parameters 
    s_s{i} = [s_s{i};s_c(2:end)];
    ps{i} = [ps{i}; ps_L(2:end,:)]; 
    Rs{i} = cat(3,Rs{i},Rs_L(:,:,2:end));
    ns{i} = [ns{i}; ns_L(2:end,:)];
    ms{i} = [ms{i}; ms_L(2:end,:)];
    
    %Calcualte n,m values at end effector for secondary rod (global)
    n_L(:,i) = ns{i}(end,:)';
    m_L(:,i) = ms{i}(end,:)';
    
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
residual = [E1(:);E2(:);E3;E4;E_inter(:);E5(:);E6(:);E7;E8];
end

