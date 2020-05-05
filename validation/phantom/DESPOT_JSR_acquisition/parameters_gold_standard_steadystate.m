%% Estimate gold standard tissue parameters
%David Leitao (david.leitao@kcl.ac.uk); 23-04-20

%Requires the handles to the SPGR and bSSFP steady-state signal functions 
%and its derivatives -> needs path to library folder

clearvars; close all; clc;

%% Load signal functions and define acquisition settings used
% vector u is the acquisition settings and vector p is the tissue parameters
%   u = [FA, TR, TE/RF] -> TE for SPGR sequences, RF for bSSFP sequences
%   p = [T1, T2, M0, P0, B0]

TEspgr = 2.3; %[ms]
Trfe   = 0.3; %[ms] RF pulse duration

% load handles to calculate steady-state signals considering finite RF
% pulses
load('hnd_SignalFunctions_FiniteRFpulses_SteadyState.mat');
Sspgr  =@(u,p) func_signal_SPGR (p(5),u(1),p(3),p(4),p(1),p(2),u(3),u(2));
Sbssfp =@(u,p) func_signal_bSSFP(p(5),u(1),p(3),p(4),u(3),p(1),p(2),u(2),Trfe);

% SPGR acquisition settings:
%              FA           TR
acq_set{1} = [[deg2rad(12), 18];                    %#1
              [deg2rad(8),  18];                    %#2 
              [deg2rad(6),  25];                    %#3
              [deg2rad(2),  35]];                   %#4

% bSSFP acquisition settings:
%              FA           TR  RF phase cycling
acq_set{2} = [[deg2rad(60), 7,  deg2rad(-70)];      %#1
              [deg2rad(60), 7,  deg2rad(70)];       %#2
              [deg2rad(10), 7,  deg2rad(-110)];     %#3
              [deg2rad(10), 7,  deg2rad(110)];      %#4
              [deg2rad(5),  7,  deg2rad(180)];      %#5
              [deg2rad(50), 7,  deg2rad(90)]];      %#6

%% Load data 

% load B1 profile for correcting nominal flip angle
load('B1_profile_gold_standard')

% load data
load('spgr_bssfp_ch2')

x_range = 1:size(imgs{1},1);
n_pts   = 1000; %no. of time-points used to averaging
            
T1_gs = zeros(numel(x_range),1);
T2_gs = zeros(numel(x_range),1);
M0_gs = zeros(numel(x_range),1);
B0_gs = zeros(numel(x_range),1);
P0_gs = zeros(numel(x_range),1);

noise_std = zeros(numel(x_range),1);

%% Fit gold standard values of [T1,T2,M0,P0,B0]

% use all data to fit the gold standard tissue parameters    
nSPGR  = 4;
nbSSFP = 6;
nSS    = nSPGR + nbSSFP;   


FA = zeros(nSS,1);
TR = zeros(nSS,1);
TE = zeros(nSPGR,1);
RF = zeros(nbSSFP,1);

for xx=1:numel(x_range)

    meas_signal_complex = zeros(nSS, n_pts);
    
    B1factor = B1_profile(x_range(xx));

    for nn=1:nSPGR
        FA(nn)           = B1factor .* acq_set{1}(nn,1);
        TR(nn)           = acq_set{1}(nn,2);
        TE(nn)           = TEspgr;
        meas_signal_complex(nn,:) = imgs{nn}(x_range(xx),end-n_pts+1:end);
        
        [C50] = CorMxySpoiling(rad2deg(FA(nn)),TR(nn));
        FA(nn) = C50 * FA(nn);
    end
    for nn=1:nbSSFP
        FA(nn+nSPGR)     = B1factor .* acq_set{2}(nn,1);
        TR(nn+nSPGR)     = acq_set{2}(nn,2);
        RF(nn)           = acq_set{2}(nn,3);
        meas_signal_complex(nn+nSPGR,:) = imgs{nn+4}(x_range(xx),end-n_pts+1:end);
    end

    meas_signal = [real(meas_signal_complex); imag(meas_signal_complex)];
    
    avg_meas_signal =  mean(meas_signal,2);

    noise_std(xx) = mean(std(meas_signal,[],2));

    u_cur = [FA; TR; TE; RF];
    f =@(p_est) SteadyStateSignals(Sspgr,Sbssfp,nSPGR,nbSSFP,u_cur,p_est);

    %    [T1    T2    M0      P0    B0]
    x0 = [1500  150   2000    0     0]; %<- starting point
    lb = [0     0     0      -pi   -Inf];
    ub = [Inf   Inf   Inf     pi    Inf];

    Fmin = @(x) (f(x) - avg_meas_signal);
    
    opts =  optimoptions('lsqnonlin',...
        'Display','none',...
        'MaxFunctionEvaluations',Inf,...
        'MaxIterations',1e4,...
        'FunctionTolerance',1e-8);
    
    x_est = lsqnonlin(Fmin,x0,lb,ub,opts);

    T1_gs(xx) = x_est(1);
    T2_gs(xx) = x_est(2);
    M0_gs(xx) = x_est(3);
    P0_gs(xx) = x_est(4);
    B0_gs(xx) = x_est(5);

    disp(['Spatial locations: ',num2str(xx),'/',num2str(numel(x_range))])

end

%% Save gold standard values

save('gold_standard_parameters','T1_gs','T2_gs','M0_gs','P0_gs','B0_gs','noise_std')

%% Plot fitted values

figure;
set(gcf,'Units','normalized','Outerposition',[0 0 1 1],'Color','w')
subplot(2,3,1)
plot(x_range,T1_gs,'Linewidth',2)
title ('T_1')
ylim([0 2000]); xlim([1 140])
ylabel('T_1 (ms)'); xlabel('x position (FE direction)')
grid minor
set(gca,'fontsize',14)

subplot(2,3,2)
plot(x_range,T2_gs,'Linewidth',2)
title ('T_2')
ylim([0 200]); xlim([1 140])
ylabel('T_2 (ms)'); xlabel('x position (FE direction)')
grid minor
set(gca,'fontsize',14)

subplot(2,3,3)
plot(x_range,M0_gs,'Linewidth',2)
title ('M_0')
ylim([0 2500]); xlim([1 140])
ylabel('M_0 (ms)'); xlabel('x position (FE direction)')
grid minor
set(gca,'fontsize',14)

subplot(2,3,4)
plot(x_range,rad2deg(P0_gs),'Linewidth',2)
title ('P_0')
ylim([-180 180]); xlim([1 140])
ylabel('P_0 (deg)'); xlabel('x position (FE direction)')
grid minor
set(gca,'fontsize',14)

load('B0_profile_gold_standard')
subplot(2,3,5)
plot(x_range,B0_gs,'Linewidth',2); 
hold on; plot(1:140,B0_profile,'Linewidth',2)
title ('B_0')
ylim([-100 100]); xlim([1 140])
ylabel('B_0 (Hz)'); xlabel('x position (FE direction)')
legend('Fitted','Dual-echo SPGR')
c = lines(2);
grid minor
set(gca,'fontsize',14)

subplot(2,3,6)
plot(x_range,B1_profile,'Linewidth',2)
title('B_1')
ylabel('B_1 factor'); xlabel('x position (FE direction)')
ylim([0.6 1.2]); xlim([1 140])
grid minor
set(gca,'fontsize',14)