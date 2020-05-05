function [ mask ] = spiral_trajectory(R, N, theta0)
%[ mask ] = spiral_trajectory(R, N, theta0)
%   Designs spiral in the cartesian grid with uniform radial density and
%   for a radial under-sampling factor of R; theta0 defines the initial 
%   phase of the spiral.
%   This design follows the concepts from Glover 1999 "Simple Analytic
%   Spiral K-space Algorithm" and Pipe and Zwart 2014 (DOI:10.1002/mrm.24675)

Nint = N/R;

theta = 0:1/(sqrt(2)*N):2*pi*Nint;
npts  = numel(theta);
kr    = linspace(0, N, npts);

k = kr .* exp(1i * (theta+theta0));

k = round(k) + (N/2+1)*(1+1i); %recentre spiral around [N/2, N/2]

% remove indexes outside the FOV
idx2rm = real(k)<1 | imag(k)<1 | real(k)>N | imag(k)>N;
k(idx2rm) = []; 

k = unique(k,'stable');
kx = real(k);
ky = imag(k)-1;
coord = kx+ky*N;

mask = zeros(N);
mask(coord) = 1;

end