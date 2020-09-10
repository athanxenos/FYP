function [v] = inv_hat(R)
%Maps SO3 matrics to R3
%Inputs:
%R - 3x3 rotation matrix

%Outputs:
%v - vector in R3

v = [R(3,2); R(1,3); R(2,1)];
end

