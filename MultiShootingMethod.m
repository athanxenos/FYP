function [residual] = MultiShootingMethod(guess)

global K_se
global K_bt
global v_ref

global d
global R_disc
global p_disc

global pb
global ps

global n
global r
global F_end
global M_end
global F_disc
global M_disc

global disc_normal
global end_normal
global pb_L


%Extract v,u sections from vector
v_guess = guess(1:30);
u_guess = guess(31:52);

%Extract backbone values at base and disc
vb0 = v_guess(1:3);
vbd = v_guess(4:6);
ub0 = u_guess(1:3);
ubd = u_guess(4:6);

%Setup secondary rod matrices
us0 = zeros(3,n);
usd = zeros(3,n);

%Extract v,u values at base and disc for secondary rods
vs0 = reshape(v_guess(7:18),[3,n]);
us0(1:2,:) = reshape(u_guess(7:14),[2,n]);
vsd = reshape(v_guess(19:30),[3,n]);
usd(1:2,:) = reshape(u_guess(15:22),[2,n]);

s_disc = guess(53:56);

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
nb_d = R_disc*K_se*(vb(end,:)'-v_ref);
mb_d = R_disc*K_bt*ub(end,:)';

%Calculate n,m guesses at first disc for central backbone (global)
nb_dguess = R_disc*K_se*(vbd-v_ref);
mb_dguess = R_disc*K_bt*ubd;

%Initialise constraint/error variables
pos_d = zeros(3,n);
ori_d = zeros(3,n);
n_d = zeros(3,n);
m_d = zeros(3,n);
md_sum = zeros(3,1);
nd_guess = zeros(3,n);
md_guess = zeros(3,n);
E_inter = zeros(1,n);
E1 = zeros(2,n);
E2 = zeros(2,n);

%Integrate secondary rods until they intersect first disc
for i=1:n
    
    %Integrate secondary rods to first disc
    [ps{i},Rs{i},vs{i},us{i},s_s{i}] = RodODE_Eval(r(:,i),Rs0,vs0(:,i),us0(:,i),0,s_disc(i));
    
    %Calculate disc intersection error
    plane = R_disc.'*(ps{i}(end,:)-p_disc)';
    E_inter(:,i) = plane(end);
    
    %Calculate positional constraint at disc (local)
    pos_d(:,i) = R_disc'*(ps{i}(end,:)-p_disc)'-r(:,i);
    
    %Positonal Error
    E1(:,i) = pos_d(1:2,i);
    
    %Calculate orientation constraint at disc
    ori_d(:,i) = inv_hat(logm(Rs{i}(:,:,end)'*R_disc));
    
    %Orientation Error
    E2(:,i) = ori_d(1:2,i);
    
    %Calculate n,m at disc (global)
    n_d(:,i) = Rs{i}(:,:,end)*K_se*(vs{i}(end,:)'-v_ref);
    m_d(:,i) = Rs{i}(:,:,end)*K_bt*us{i}(end,:)';
    
    %Calculate n,m guesses at first disc (global)
    nd_guess(:,i) = Rs{i}(:,:,end)*K_se*(vsd(:,i)-v_ref);
    md_guess(:,i) = Rs{i}(:,:,end)*K_bt*usd(:,i);

    %Sum moments at disc
    md_sum = md_sum + cross(ps{i}(end,:)',(nd_guess(:,i) - n_d(:,i))) + md_guess(:,i) - m_d(:,i);
end

%Sum forces at disc end
nd_minus = sum(n_d,2) + nb_d;

%Add guesses for forces at disc start
nd_plus = sum(nd_guess,2) + nb_dguess;

%Force Equilibrium Error
E3 = nd_plus + nd_minus - F_disc;

%Moment Equilibrium Error
E4 = md_sum + cross(p_disc',(nb_dguess - nb_d - F_disc)) + mb_dguess - mb_d - M_disc;

%Integrate central backbone from first disc to end effector
[pb_L,Rb_L,vb_L,ub_L,s_L] = RodODE_Eval(pb(end,:)',Rb(:,:,end),vbd,ubd,d(1),d(2));

%Concatenate s,p,R,v,u parameters 
s = [s;s_L(2:end)];
pb = [pb;pb_L(2:end,:)];
Rb = cat(3,Rb,Rb_L(:,:,2:end));
vb = [vb;vb_L(2:end,:)];
ub = [ub;ub_L(2:end,:)];

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
    [ps_L,Rs_L,vs_L,us_L,s_c] = RodODE_Eval(ps{i}(end,:)',Rs{i}(:,:,end),vsd(:,i),usd(:,i),s_disc(i),d(2)); 
    
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
    E5(:,i) = Rb_L'*(ps{i}(end,:)' - pb_L)-r(:,i);
    
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

