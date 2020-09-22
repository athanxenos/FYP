function [PosDisc,OriDisc,FDisc,MDisc,PosEnd,OriEnd,FEnd,MEnd] = residual_extract(residual)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
PosDisc = residual(1:8);
OriDisc = residual(9:16);
FDisc = residual(17:19);
MDisc = residual(20:22);
PosEnd = residual(23:34);
OriEnd = residual(35:46);
FEnd = residual(47:49);
MEnd = residual(50:52);
end

