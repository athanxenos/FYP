function [u_back_base, u_back_disc, u_sec_base, u_sec_disc] = guess_extract(guess)
%Extracts guess into u components for backbone and secondayr rods

n=4;
u = guess(23:52);
u_back_base = u(1:3);
u_back_disc = u(4:6);
u_sec_base = reshape(u(7:18),[3,n]);
u_sec_disc = reshape(u(19:end),[3,n]);
end

