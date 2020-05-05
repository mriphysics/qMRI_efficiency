%% Estimation of B1 using AFI method
%David Leitao (david.leitao@kcl.ac.uk); 23-04-20

clearvars; close all; clc;

%% Load AFI data and estimate gold standard B1+ values

load('AFI_SPGR')

echo1 = imgs(:,1:2:end);
echo2 = imgs(:,2:2:end);

a_nominal = deg2rad(60);
TR1 = 25;
TR2 = 100;
n = TR2/TR1;

r = double(abs(echo2./echo1));

B1_est = real(acos((r.*n-1)./(n-r))./a_nominal);

figure; 
imagesc(B1_est)
caxis([0.8 1.7]); colorbar;

% average final estimate over the last 800 time-points
B1_profile = mean(B1_est(:,end-800:end),2);
plot(B1_profile)

save('B1_profile_gold_standard','B1_profile')
