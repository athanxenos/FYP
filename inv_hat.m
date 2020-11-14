function [v] = inv_hat(R)
%inv_hat
%Written by Athan Xenos

%Maps SO3 matrices to R3 vector
%Inputs:
%R - 3x3 rotation matrix

%Outputs:
%v - vector in R3

v = [R(3,2); R(1,3); R(2,1)];
end

