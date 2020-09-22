function [u_back_base, u_back_disc, u_sec_base, u_sec_disc] = guess_extract(guess)
%Extracts guess into u components for backbone and secondayr rods

n=4;
u_sec_base = zeros(3,n);
u_sec_disc = zeros(3,n);

u = guess(31:end);
u_back_base = u(1:3);
u_back_disc = u(4:6);
u_sec_base(1:2,:) = reshape(u_guess(7:14),[2,n]);
u_sec_disc(1:2,:) = reshape(u_guess(15:22),[2,n]);
end

