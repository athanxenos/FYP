function [u_hat] = hat(u)
%hat
%Written by Athan Xenos

%Hat operator function
%Converts R3 vector to 3x3 matrix, R6 vector to 4x4 matrix

n=length(u);
if n==3
    
    ux = u(1);
    uy = u(2);
    uz = u(3);
    u_hat = [0 -uz uy; uz 0 -ux; -uy ux 0];
    
elseif n==6
    
    vx = u(1);
    vy = u(2);
    vz = u(3);
    ux = u(4);
    uy = u(5);
    uz = u(6);
    u_hat = [0 -uz uy vx; uz 0 -ux vy; -uy ux 0 vz; 0 0 0 0];
    
else 
    error('Incorrect input argument dimension (must be length 3 or 6)')
end

