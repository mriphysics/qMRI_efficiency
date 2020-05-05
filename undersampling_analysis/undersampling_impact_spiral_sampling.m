%% Monte Carlo simulation of the signals with spiral sampling
%David Leitao (david.leitao@kcl.ac.uk); 23-04-20

clearvars; close all; clc;

%% Run Monte Carlo simulation for several under-sampling factors

N = 128;
nimgs = 1e5;

addThermalNoise = 1; %<- change this flag to add or not thermal noise
SNR = 30/N; %<-change this value to try different SNR levels; check line 10

% use Shepp-Logan phantom
ytrue = phantom(N);
Ytrue = fftshift(fft2(ytrue));

ymeas  = zeros(N,N,nimgs);

% under-sampling factors
R = [1 2 4 6 8 10 12 16 20 24 30 36 42 48];
rng('default')

% rotate the spiral with the tiny gold angle (N=7), following Wundrak et 
% al. 2016 (DOI:10.1002/mrm.25831)
TGA7 = pi / ( (1+sqrt(5))/2 + 7 -1 );

load('spiral_DCF.mat')

nn = 1;
for rr=R
    for ii=1:nimgs
        if addThermalNoise
            thermal_noise = randn(N)/SNR;
            Ytrue0 = Ytrue + thermal_noise; 
        else
            Ytrue0 = Ytrue;
        end
        
        % design spiral in cartesian space to define measured k-space
        phi = (ii-1) * TGA7;
        mask = spiral_trajectory(rr, N, phi);
        
        Ymeas = Ytrue0 .* mask .* DCF(:,:,nn); %applying k-space sampling aand DCF

        %reconstruction of zero-filled image
        ymeas(:,:,ii) = ifft2(ifftshift(Ymeas));

    end
    
    aliasing.avg(:,:,nn) = mean(ymeas,3);
    aliasing.std(:,:,nn) = std(ymeas,[],3);
    
    nn = nn + 1;
end

%% Save results

save('spiral_sampling_SNR30',...
    'aliasing','nimgs','N','R','SNR')

% save('spiral_sampling_SNR50',...
%     'aliasing','nimgs','N','R','SNR')

% save('spiral_sampling_SNR80',...
%     'aliasing','nimgs','N','R','SNR')

% save('spiral_sampling_SaR',...
%     'aliasing','nimgs','N','R')



