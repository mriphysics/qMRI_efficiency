%% Calculation of the Density Compensation Function for spiral sampling
%David Leitao (david.leitao@kcl.ac.uk); 23-04-20

% DCF is calculated by generating many spiral for each under-sampling
% factor to get the sampling density and then matching it to a uniform
% density

clearvars; close all; clc;

%%

N = 128;
nimgs = 1e5;

% under-sampling factors
R = [1 2 4 6 8 10 12 16 20 24 30 36 42 48];
rng('default')

Ysamp  = zeros(N, N, numel(R)); %<- records the sampled points of every trial

% rotate the spiral with the tiny gold angle (N=7), following Wundrak et 
% al. 2016 (DOI:10.1002/mrm.25831)
TGA7 = pi / ( (1+sqrt(5))/2 + 7 -1 );

nn = 1;
for rr=R
    for ii=1:nimgs  
        phi = (ii-1) * TGA7; 
        mask = spiral_trajectory(rr, N, phi);

        Ysamp(:,:,nn) = Ysamp(:,:,nn) + mask;
    end
    nn = nn + 1;
end

DCF = nimgs./Ysamp; % overall density of each k-space location
for nn=1:numel(R)
    %DCF scales each k-space location to be sampled with density 1/R
    DCF(:,:,nn) = DCF(:,:,nn) / R(nn); 
end

save('spiral_DCF','DCF','R')
