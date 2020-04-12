function [xArc, yArc] = plotArc(x, r)
%Returns (x,y) co-ordinates for arc circle centred at (r,0) with radius r

xArc = linspace(0,x);
yArc = sqrt(r^2-(xArc-r).^2);
end

