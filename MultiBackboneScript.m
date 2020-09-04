%Cosserat Model Script based on 2017 paper
%Model multiple backbones
clear all
close all
clc 

%Define global variables for model
global K_se
global K_bt
global v_ref

%Rod Parameters
L = 0.25; %Arclength of all rods (m)
rad = 0.0005; %Radius of all rods (m)(approx 0.5mm)
Area = pi*rad^2; %Area of cross section (m^2)

%Material parameters
EY = 210*10^9; %Young's Modulus (Pa)
GY = 80*10^9;  %Shear Modulus (Pa) G=E/2(1+v) (v=0.3125)

%Second moments of area of cross section
Ixx = pi*rad^4/4;
Iyy = pi*rad^4/4;
Izz = Ixx + Iyy; %Polar moment of inertia

%Stiffness Matrices
K_se = diag([GY*Area GY*Area EY*Area]); %Shear/extension stiffness matrix
K_bt = diag([EY*Ixx EY*Iyy EY*Izz]); %Bending/torsion stiffness matrix

%Rod Parameters
rad_s = 0.015; %Radial location of secondary rods from central backbone (m)(approx 10mm)
n_s = 4; %Number of secondary backbones
ps0 = [0 -rad_s 0 rad_s; rad_s 0 -rad_s 0; 0 0 0 0]; %Radial coordinate profile of secondary backbones

%Disc Parameters
n_d = 2;  %Number of discs after base(including end effector)
d = [L/2;L]; %Disc locations on central backbone

%Reference Parameters
%Linear/angular rate of change of frame in reference state
v_ref = [0;0;1];

%/////////// Model Variables ////////////
%Guess initial conditions for v,u of central backbone
v_guess=[0;0;1]; %Linear rate of change of frame
u_guess=[-1;0;0]; %Angular rate of change of frame
%////////////////////////////////////////

%Define initial conditions for central backbone
pb0 = [0;0;0];
Rb0 = eye(3);
vb0 = v_guess;
ub0 = u_guess;

%Integrate central backbone upto first disc
[pb,Rb,vb,ub,s] = RodIntegrate(pb0,Rb0,vb0,ub0,0,d(1));

%Extract normal to disc as z axis of backbone R matrix
R_disc = Rb(:,:,end);
p_disc = pb(end,:);
disc_normal = R_disc(:,3);

%Initialise variables for secondary backbones
Rs0 = zeros(3,3,n_s);
vs0 = zeros(3,n_s);
us0 = zeros(3,n_s);

ps = cell(1,n_s);
us = cell(1,n_s);
vs = cell(1,n_s);
Rs = cell(1,n_s);
s_s = cell(1,n_s);

%Step size to look before and after disc for rod intersection (20mm)
step = 0.015;

for i=1:n_s
    
    %Define initial conditions at base
    Rs0(:,:,i) = eye(3);
    vs0(:,i) = v_guess;
    us0(:,i) = u_guess;
    
    %Integrate secondary rods to just before first disc (Section A)
    [ps_a,Rs_a,vs_a,us_a,s_a] = RodIntegrate(ps0(:,i),Rs0(:,:,i),vs0(:,i),us0(:,i),0,d(1)-step);
    
    %Integrate secondary rods from step before first disc to step after
    %first disc
    [ps_b,Rs_b,vs_b,us_b,s_b] = RodIntegrate(ps_a(end,:)',Rs_a(:,:,end),vs_a(end,:)',us_a(end,:)',d(1)-step, d(1)+step);
    
    %Combine s and p parameters for both sections
    s_s{i} = [s_a; s_b(2:end)];
    ps{i} = [ps_a; ps_b(2:end,:)]; 
    Rs{i} = cat(3,Rs_a,Rs_b(:,:,2:end));
    vs{i} = [vs_a; vs_b(2:end,:)];
    us{i} = [us_a; us_b(2:end,:)];
end

%Combine section lengths
n = length(s_s{1});
n_a = length(s_a);

%Find index of first disc location
disc_index = find(s_s{i}==L/2);
disc_check = ones(n,n_s);

%Loop through each secondary backbone
for i=1:n_s
    for j=n_a:n     %Loop through each entry in section B
        plane_vec = R_disc.'*(ps{i}(j,:)-p_disc)';  %Check intersection condition
        disc_check(j,i) = plane_vec(end);
         if norm(disc_check(j,i))<0.0005 && j==disc_index %Special case for intersection of disc at central backbone location
             disc_check(j,i)=0;
             fprintf('Intersection found at central backbone location\n');
             break  %Break once intersection found
         elseif norm(disc_check(j,i))<0.0005    %Search for intersection of disc within tolerance
             disc_check(j,i)=0;
             fprintf('Intersection found at disc\n');
             break  %Break once intersection found
         end
        
     end
end

%Find locations of intersection points
[inter, col] = find(disc_check == 0);

%Shorten arclength to intersection point
%s_s = s_s(1:max(s_inter),:);
if length(inter) < 4
    fprintf('Error: 4 disc intersection points not found');
end

%Concatenate secondary rods upto first disc
for i=1:n_s
    ps{i} = ps{i}(1:inter(i),:);
    vs{i} = vs{i}(1:inter(i),:);
    us{i} = us{i}(1:inter(i),:);
    Rs{i} = Rs{i}(:,:,1:inter(i));
    s_s{i} = s_s{i}(1:inter(i));
    
    mid(i) = s_s{i}(inter(i));
end

%Integrate central backbone from first disc to end effector
[pb_end,Rb_end,vb_end,ub_end,s_end] = RodIntegrate(pb(end,:)',Rb(:,:,end),vb(end,:)',ub(end,:)',d(1),d(2));
%[pb,Rb,vb,ub,s] = RodIntegrate(pb0,Rb0,vb0,ub0,0,d(2));
s = [s;s_end(2:end)];
pb = [pb;pb_end(2:end,:)];
Rb = cat(3,Rb,Rb_end(:,:,2:end));
vb = [vb;vb_end(2:end,:)];
ub = [ub;ub_end(2:end,:)];

%Integrate from first disc to end effector
for i=1:n_s
    
    %Integrate secondary rods from first disc to end effector (Section C)
    [ps_c,Rs_c,vs_c,us_c,s_c] = RodIntegrate(ps{i}(end,:)',Rs{i}(:,:,end),vs{i}(end,:)',us{i}(end,:)',mid(i),d(2)); 
    
    %Combine s,p,R,v,u parameters for all sections
    s_s{i} = [s_s{i};s_c(2:end)];
    ps{i} = [ps{i}; ps_c(2:end,:)]; 
    Rs{i} = cat(3,Rs{i},Rs_c(:,:,2:end));
    vs{i} = [vs{i}; vs_c(2:end,:)];
    us{i} = [us{i}; us_c(2:end,:)];
end


%Calculate arclength to check solution feasibility
arc = arclength(pb(:,1),pb(:,2),pb(:,3));

centre_index = find(s==L/2);
d_centre = pb(centre_index,:);
end_effector = pb(end,:);
end_normal = Rb(:,3,end);

%Plot solution for central backbone
hold on

plot3(pb(:,1),pb(:,2),pb(:,3),'b');
%plot3(pb(centre_index,1),pb(centre_index,2),pb(centre_index,3),'ok');


%Loop and plot secondary backbones
for i = 1:n_s
    plot3(ps{i}(:,1),ps{i}(:,2),ps{i}(:,3),'r');
    %plot3(ps{i}(inter(i),1),ps{i}(inter(i),2),ps{i}(inter(i),3),'ok');
end

plotCircle3D(d_centre,disc_normal',rad_s)
plotCircle3D(end_effector,end_normal',rad_s)

xlabel('x');
ylabel('y');
zlabel('z');
grid on
axis([-L,L,-L,L,-L,L]);
title(['Arclength is ',num2str(arc)])