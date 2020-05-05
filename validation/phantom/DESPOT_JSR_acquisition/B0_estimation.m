%% Estimation of B0 using a dual echo SPGR
%David Leitao (david.leitao@kcl.ac.uk); 23-04-20

clearvars; close all; clc;

%% Load Dual echo SPGR data and calculate gold standard B0 values

load('DualEchoSPGR')

echo1 = imgs(:,1:2:end);
echo2 = imgs(:,2:2:end);

phase1 = angle(echo1);
phase2 = angle(echo2);

figure;
subplot(2,1,1)
imagesc(phase1)
subplot(2,1,2)
imagesc(phase2)

dTE = 2.3e-3; %echo time difference between the two signals

B0_est = angle(echo2./echo1)/(2*pi*dTE);

figure; 
imagesc(B0_est)

% average final estimate over the last 2000 time-points
B0_profile = mean(B0_est(:,end-2000:end),2); 
plot(B0_profile)

save('B0_profile_gold_standard','B0_profile')
