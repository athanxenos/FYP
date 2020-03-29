function [T] = TMatrix(s,k,phi)
%Homogenous Transformation Matrix for Jones Kinematic Model
%Maps configuration space to task space
%Phi in radians

T = [cos(phi) -sin(phi)*cos(s*k) sin(phi)*sin(s*k) sin(phi)*(1-cos(s*k))/k;
    sin(phi) cos(phi)*cos(s*k) -cos(phi)*sin(s*k) -cos(phi)*(1-cos(s*k))/k;
    0 sin(s*k) cos(s*k) sin(s*k)/k;
    0 0 0 1];  
end

