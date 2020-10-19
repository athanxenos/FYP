function [residual] = MultiShootingMethod(guess)
%Function that evaluates Multi Backbone ODE's for input guess and returns
%residual of position/orientation/equilibrium conditions

%Inputs:
% guess - 64x1 vector with initial guess for n,m,s_disc (52 length minimum)
% All v/u values are initialised as [0;0;1] and [0;0;0]
% s_disc intersections are initialised for straight rod position

%Outputs:
% residual - returns error vector for given guess (56 length minimum)

%Define global variables for model
%Model Parameters
global d
global n
global nd
global n_mid
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

residual = zeros(30+26*nd,1);

%Define initial conditions for central backbone
pb0 = [0;0;0];
Rb0 = eye(3);

%Initialise variables for secondary backbones
Rs0 = repmat(eye(3),1,1,4);
ps0 = r;

%Initialise cells for secondary backbones
ps = cell(1,n);
ns = cell(1,n);
ms = cell(1,n);
Rs = cell(1,n);
s_s = cell(1,n);

disc_normal = zeros(nd,3);
R_disc = zeros(3,3,nd+1);
R_disc(:,:,1) = Rb0;
p_disc = zeros(nd,3);
p_disc(1,:) = pb0;
s_disc = zeros(nd+1,4);

for j =1:nd
    
    %Extract guess for central rod
    nb0 = guess(1+30*(j-1):3+30*(j-1));
    mb0 = guess(4+30*(j-1):6+30*(j-1));
    
    nb_Dguess = guess(1+30*j:3+30*j);
    mb_Dguess = guess(4+30*j:6+30*j);
    
    
    
    %Integrate central backbone upto first disc
    [pb_temp,Rb_temp,nb_temp,mb_temp,s_temp] = RodODE_Eval_force(p_disc(j,:)',R_disc(:,:,j),nb0,mb0,d(j),d(j+1));
    
    if j==1
        s = s_temp;
        pb = pb_temp;
        Rb = Rb_temp;
        nb = nb_temp;
        mb = mb_temp;
    else
   
        %Concatenate s,p,R,v,u parameters 
        s = [s;s_temp(2:end)];
        pb = [pb;pb_temp(2:end,:)];
        Rb = cat(3,Rb,Rb_temp(:,:,2:end));
        nb = [nb;nb_temp(2:end,:)];
        mb = [mb;mb_temp(2:end,:)];
    end
    
    %Extract normal to disc as z axis of backbone R matrix
    R_disc(:,:,j+1) = Rb(:,:,end);
    p_disc(j+1,:) = pb(end,:);
    disc_normal(j,:) = R_disc(:,3,j+1);

    %Calculate n,m at disc for central backbone (global)
    nb_D = nb(end,:)';
    mb_D = mb(end,:)';

    %Initialise constraint/error variables
    pos_D = zeros(3,n);
    ori_D = zeros(3,n);
    n_D = zeros(3,n);
    m_D = zeros(3,n);
    mD_sum = zeros(3,1);
    E1 = zeros(2,n);
    E2 = zeros(2,n);
    E_inter = zeros(1,n);
    ns_Dguess = zeros(3,n);
    ms_Dguess = zeros(3,n);
    
    s_disc(j+1,:) = guess(30*(nd+1)+1+4*(j-1):30*(nd+1)+4+4*(j-1));
    
    %Integrate secondary rods until they intersect first disc
    for i=1:n
        
        ns0 = guess(7+30*(j-1)+3*(i-1):9+30*(j-1)+3*(i-1));
        ms0 = guess(19+30*(j-1)+3*(i-1):21+30*(j-1)+3*(i-1));
        
        ns_Dguess(:,i) = guess(37+30*(j-1)+3*(i-1):39+30*(j-1)+3*(i-1));
        ms_Dguess(:,i) = guess(49+30*(j-1)+3*(i-1):51+30*(j-1)+3*(i-1));
        
        
        %Integrate secondary rods to first disc
        [ps_temp,Rs_temp,ns_temp,ms_temp,s_temp] = RodODE_Eval_force(ps0(:,i),Rs0(:,:,i),ns0,ms0,s_disc(j,i),s_disc(j+1,i));
        
        %Concatenate s,p,R,v,u parameters 
        s_s{i} = [s_s{i};s_temp(2:end)];
        ps{i} = [ps{i}; ps_temp(2:end,:)]; 
        Rs{i} = cat(3,Rs{i},Rs_temp(:,:,2:end));
        ns{i} = [ns{i}; ns_temp(2:end,:)];
        ms{i} = [ms{i}; ms_temp(2:end,:)];
        
        Rs0(:,:,i) = Rs{i}(:,:,end);
        ps0(:,i) = ps{i}(end,:)';
        
        %Calculate disc intersection error
        intersect = R_disc(:,:,j+1).'*(ps{i}(end,:)-p_disc(j+1,:))';
        E_inter(:,i) = intersect(end);

        %Calculate positional constraint at disc (local)
        pos_D(:,i) = R_disc(:,:,j+1)'*(ps{i}(end,:)-p_disc(j+1,:))' - r(:,i);

        %Positonal Error
        E1(:,i) = pos_D(1:2,i);

        %Calculate orientation constraint at disc
        ori_D(:,i) = inv_hat(logm(Rs{i}(:,:,end)'*R_disc(:,:,j+1)));

        %Orientation Error
        E2(:,i) = ori_D(1:2,i);

        %Calculate n,m at disc (global)
        n_D(:,i) = ns{i}(end,:)';
        m_D(:,i) = ms{i}(end,:)';

        %Sum moments at disc
        mD_sum = mD_sum + cross(ps{i}(end,:)',(-ns_Dguess(:,i) + n_D(:,i))) - ms_Dguess(:,i) + m_D(:,i);
    end

    %Sum forces exiting disc
    nd_minus = sum(n_D,2) + nb_D;

    %Sum forces entering disc
    nd_plus = sum(ns_Dguess,2) + nb_Dguess;

    if j == n_mid
        
        %Force Equilibrium Error
        E3 = -nd_plus + nd_minus - F_disc;

        %Moment Equilibrium Error
        E4 = mD_sum + cross(p_disc(j+1,:)',(-nb_Dguess + nb_D - F_disc)) - mb_Dguess + mb_D - M_disc;
        
    else
        %Force Equilibrium Error
        E3 = -nd_plus + nd_minus;

        %Moment Equilibrium Error
        E4 = mD_sum + cross(p_disc(j+1,:)',(-nb_Dguess + nb_D)) - mb_Dguess + mb_D;
    end
    
    residual(1+26*(j-1):26+26*(j-1)) = [E1(:);E2(:);E3;E4;E_inter(:)];
    
end

%% //////////// Last Disc to End Effector ////////////
%Integrate central backbone from first disc to end effector
[pb_end,Rb_end,nb_end,mb_end,s_L] = RodODE_Eval_force(pb(end,:)',Rb(:,:,end),nb_Dguess,mb_Dguess,d(end-1),d(end));

%Concatenate s,p,R,v,u parameters 
s = [s;s_L(2:end)];
pb = [pb;pb_end(2:end,:)];
Rb = cat(3,Rb,Rb_end(:,:,2:end));
nb = [nb;nb_end(2:end,:)];
mb = [mb;mb_end(2:end,:)];

%Find end effector position/orientation
pb_L = pb(end,:)';       
Rb_L =  Rb(:,:,end);
end_normal = Rb_L(:,3);

%Calculate n,m at end effector for central backbone (global)
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
    
    ns_Dguess = guess(7+30*nd+3*(i-1):9+30*nd+3*(i-1));
    ms_Dguess = guess(19+30*nd+3*(i-1):21+30*nd+3*(i-1));
    
    %Integrate secondary rods from first disc to end effector 
    [ps_L,Rs_L,ns_L,ms_L,s_c] = RodODE_Eval_force(ps{i}(end,:)',Rs{i}(:,:,end),ns_Dguess,ms_Dguess,s_disc(end,i),d(end)); 
    
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
residual(1+26*nd:30+26*nd) = [E5(:);E6(:);E7;E8];

end

