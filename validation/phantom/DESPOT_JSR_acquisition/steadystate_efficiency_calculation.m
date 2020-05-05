%% Calculation and comparison of experimental and theoretical efficiency
%David Leitao (david.leitao@kcl.ac.uk); 23-04-20

%Requires the handles to the SPGR and bSSFP steady-state signal functions 
%and its derivatives -> needs path to library folder

clearvars; close all; clc;

%% Load signal functions and its derivatives;
% vector u is the acquisition settings and vector p is the tissue parameters
%   u = [FA, TR, TE/RF] -> TE for SPGR sequences, RF for bSSFP sequences
%   p = [T1, T2, M0, P0, B0]

TEspgr = 2.3; %[ms]
Trfe   = 0.3; %[ms]

% load handles to calculate steady-state signals considering finite RF
% pulses
load('hnd_SignalFunctions_FiniteRFpulses_SteadyState.mat');
Sspgr  =@(u,p) func_signal_SPGR (p(5),u(1),p(3),p(4),p(1),p(2),u(3),u(2));
Sbssfp =@(u,p) func_signal_bSSFP(p(5),u(1),p(3),p(4),u(3),p(1),p(2),u(2),Trfe);

% load handles to calculate signal derivatives considering finite RF pulses
load('hnd_SignalDerivatives_FiniteRFpulses_SteadyState.mat')
dmdT1_SPGR  =@(u,p) func_dmdT1_SPGR (p(5),u(1),p(3),p(4),p(1),p(2),u(3),u(2)); 
dmdT1_bSSFP =@(u,p) func_dmdT1_bSSFP(p(5),u(1),p(3),p(4),u(3),p(1),p(2),u(2),Trfe);
dmdT2_SPGR  =@(u,p) func_dmdT2_SPGR (p(5),u(1),p(3),p(4),p(1),p(2),u(3),u(2));
dmdT2_bSSFP =@(u,p) func_dmdT2_bSSFP(p(5),u(1),p(3),p(4),u(3),p(1),p(2),u(2),Trfe);
dmdM0_SPGR  =@(u,p) func_dmdM0_SPGR (p(5),u(1),p(4),p(1),p(2),u(3),u(2));
dmdM0_bSSFP =@(u,p) func_dmdM0_bSSFP(p(5),u(1),p(4),u(3),p(1),p(2),u(2),Trfe);
dmdP0_SPGR  =@(u,p) func_dmdP0_SPGR (p(5),u(1),p(3),p(4),p(1),p(2),u(3),u(2));
dmdP0_bSSFP =@(u,p) func_dmdP0_bSSFP(p(5),u(1),p(3),p(4),u(3),p(1),p(2),u(2),Trfe);
dmdB0_SPGR  =@(u,p) func_dmdB0_SPGR (p(5),u(1),p(3),p(4),p(1),p(2),u(3),u(2));
dmdB0_bSSFP =@(u,p) func_dmdB0_bSSFP(p(5),u(1),p(3),p(4),u(3),p(1),p(2),u(2),Trfe);

%% Define acquisition settings used

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


% combinations used for calculating the experimental efficiency
A5 = combnk(1:10,5);
idx5 = A5;
idx5(A5<5) = 0; idx5(A5>=5) = 1; %idx contains 1 where is a bSSFP signal
A5((sum(idx5,2)<3),:) = []; %exclude combinations with less than 2 bSSFP 
A6 = combnk(1:10,6);
idx6 = A6;
idx6(A6<5) = 0; idx6(A6>=5) = 1;
idx6 = (A6.*idx6) > 0;
A6((sum(idx6,2)<3),:) = []; %exclude combinations with less than 2 bSSFP
A7 = combnk(1:10,7);
idx7 = A7;
idx7(A7<5) = 0; idx7(A7>=5) = 1; 
idx7 = (A7.*idx7) > 0;
A7((sum(idx7,2)<3),:) = []; %exclude combinations with less than 2 bSSFP
A8 = combnk(1:10,8);
idx8 = A8;
idx8(A8<5) = 0; idx8(A8>=5) = 1; 
idx8 = (A8.*idx8) > 0;
A8((sum(idx8,2)<3),:) = []; %exclude combinations with less than 2 bSSFP and with just bSSFP
A9 = combnk(1:10,9);
idx9 = A9;
idx9(A9<5) = 0; idx9(A9>=5) = 1; 
idx9 = (A9.*idx9) > 0;
A9((sum(idx9,2)<3),:) = []; %exclude combinations with less than 2 bSSFP and with just bSSFP
A10 = 1:10;

ncomb = size(A5,1)+size(A6,1)+size(A7,1)+size(A8,1)+size(A9,1)+size(A10,1);

combinations = cell(ncomb,2);
for ii=1:ncomb
    if ii<=size(A5,1)
        aux_ii = ii;
        combinations{ii,1} = A5(aux_ii,A5(aux_ii,:)<5);
        combinations{ii,2} = A5(aux_ii,A5(aux_ii,:)>=5);
    elseif ii<=size(A5,1)+size(A6,1) && ii>size(A5,1)
        aux_ii = ii - size(A5,1);
        combinations{ii,1} = A6(aux_ii,A6(aux_ii,:)<5);
        combinations{ii,2} = A6(aux_ii,A6(aux_ii,:)>=5);
    elseif ii<=size(A5,1)+size(A6,1)+size(A7,1) && ii>size(A5,1)+size(A6,1)
        aux_ii = ii - (size(A5,1)+size(A6,1));
        combinations{ii,1} = A7(aux_ii,A7(aux_ii,:)<5);
        combinations{ii,2} = A7(aux_ii,A7(aux_ii,:)>=5);
    elseif ii<=size(A5,1)+size(A6,1)+size(A7,1)+size(A8,1) && ii>size(A5,1)+size(A6,1)+size(A7,1)
        aux_ii = ii - (size(A5,1)+size(A6,1)+size(A7,1));
        combinations{ii,1} = A8(aux_ii,A8(aux_ii,:)<5);
        combinations{ii,2} = A8(aux_ii,A8(aux_ii,:)>=5);
    elseif ii<=size(A5,1)+size(A6,1)+size(A7,1)+size(A8,1)+size(A9,1) && ii>size(A5,1)+size(A6,1)+size(A7,1)+size(A8,1)
        aux_ii = ii - (size(A5,1)+size(A6,1)+size(A7,1)+size(A8,1));
        combinations{ii,1} = A9(aux_ii,A9(aux_ii,:)<5);
        combinations{ii,2} = A9(aux_ii,A9(aux_ii,:)>=5);
    else
        aux_ii = ii - (size(A5,1)+size(A6,1)+size(A7,1)+size(A8,1)+size(A9,1));
        combinations{ii,1} = A10(aux_ii,A10(aux_ii,:)<5);
        combinations{ii,2} = A10(aux_ii,A10(aux_ii,:)>=5);
    end
end
            
%% Load data

xpnt = 60; %select voxel to calculate the efficiency
npts = 1000;

% load B1 profile for correcting nominal flip angle
load('B1_profile_gold_standard')

% load data
load('spgr_bssfp_ch2')

% load gold standard parameters 
load('gold_standard_parameters.mat')
T1_gs     = T1_gs(xpnt);
T2_gs     = T2_gs(xpnt);
M0_gs     = M0_gs(xpnt);
P0_gs     = P0_gs(xpnt);
B0_gs     = B0_gs(xpnt);
noise_std = noise_std(xpnt);

p_gs = [T1_gs T2_gs 1 P0_gs B0_gs];

experimental_eff = zeros(ncomb,2);           
theoretical_eff  = zeros(ncomb,2);           

T1_est = zeros(ncomb, npts);
T2_est = zeros(ncomb, npts);
M0_est = zeros(ncomb, npts);
B0_est = zeros(ncomb, npts);
P0_est = zeros(ncomb, npts);

%% Calculate experimental efficiency

% /!\ This section makes use of parallel acceleration: to disable it 
% comment lines 149-153 and replace parfor loop in line 196 by a for loop

if isempty(gcp('nocreate'))
    c = parcluster('local');
    c.NumWorkers = 4;
    parpool(c, c.NumWorkers);
end

%    [T1    T2    M0      P0    B0]
x0 = [1500  150   2000    0     0]; %<- starting point
lb = [0     0     0      -pi   -Inf];
ub = [Inf   Inf   Inf     pi    Inf];

for ii=1:ncomb
    nSPGR  = numel(combinations{ii,1});
    nbSSFP = numel(combinations{ii,2});
    nSS    = nSPGR + nbSSFP;   
    
    FA = zeros(nSS,1);
    TR = zeros(nSS,1);
    TE = zeros(nSPGR,1);
    RF = zeros(nbSSFP,1);

    meas_signal_complex = zeros(nSS, npts);
    
    B1factor = B1_profile(xpnt);

    for nn=1:nSPGR
        FA(nn)           = B1factor .* acq_set{1}(combinations{ii,1}(nn),1);
        TR(nn)           = acq_set{1}(combinations{ii,1}(nn),2);
        TE(nn)           = TEspgr;
        meas_signal_complex(nn,:)      = imgs{combinations{ii,1}(nn)}(xpnt,end-npts+1:end);
        
        [C50] = CorMxySpoiling(rad2deg(FA(nn)),TR(nn));
        FA(nn) = C50 * FA(nn);
    end
    for nn=1:nbSSFP
        FA(nn+nSPGR)     = B1factor .* acq_set{2}(combinations{ii,2}(nn)-4,1);
        TR(nn+nSPGR)     = acq_set{2}(combinations{ii,2}(nn)-4,2);
        RF(nn)           = acq_set{2}(combinations{ii,2}(nn)-4,3);
        meas_signal_complex(nn+nSPGR,:) = imgs{combinations{ii,2}(nn)}(xpnt,end-npts+1:end);
    end
 
    meas_signal = [real(meas_signal_complex); imag(meas_signal_complex)];

    u_cur = [FA; TR; TE; RF];
    f =@(p_est) SteadyStateSignals(Sspgr,Sbssfp,nSPGR,nbSSFP,u_cur,p_est);
    
    tic
    parfor nn=1:npts
        single_meas_signal = meas_signal(:,nn);
        Fmin = @(x) (f(x) - single_meas_signal);

        opts =  optimoptions('lsqnonlin',...
            'Display','none',...
            'MaxFunctionEvaluations',Inf,...
            'MaxIterations',1e4,...
            'FunctionTolerance',1e-8);
        [x_est,res_norm] = lsqnonlin(Fmin,x0,lb,ub,opts);
        
        T1_est(ii,nn) = x_est(1);
        T2_est(ii,nn) = x_est(2);
        M0_est(ii,nn) = x_est(3);
        P0_est(ii,nn) = x_est(4);
        B0_est(ii,nn) = x_est(5);
    end

    T1_avg = mean(mean(T1_est(ii,:))); T1_std = mean(std( T1_est(ii,:)));
    T2_avg = mean(mean(T2_est(ii,:))); T2_std = mean(std( T2_est(ii,:)));
    M0_avg = mean(mean(M0_est(ii,:))); M0_std = mean(std( M0_est(ii,:)));
    P0_avg = mean(mean(P0_est(ii,:))); P0_std = mean(std( P0_est(ii,:)));
    B0_avg = mean(mean(B0_est(ii,:))); B0_std = mean(std( B0_est(ii,:)));  
    
    experimental_eff(ii,1) = T1_gs * noise_std / ...
       (T1_std * M0_gs * sqrt(sum(TR)*1e-3));
   
    experimental_eff(ii,2) = T2_gs * noise_std / ...
       (T2_std * M0_gs * sqrt(sum(TR)*1e-3));
   
   [theoretical_eff(ii,:), ~] = analysis_acq_set_DESPOT_JSR(u_cur, p_gs, nbSSFP, ...
       dmdT1_SPGR,  dmdT2_SPGR,  dmdM0_SPGR,  dmdP0_SPGR,  dmdB0_SPGR, ...
       dmdT1_bSSFP, dmdT2_bSSFP, dmdM0_bSSFP, dmdP0_bSSFP, dmdB0_bSSFP,...
       Sspgr, Sbssfp);
   
   disp(' ')
    disp([num2str(ii),'/',num2str(ncomb)])
   disp(['Theoretical:  ',char(414),'(T1)=',num2str(theoretical_eff(ii,1),'%.3f'),'   ',char(414),'(T2)=',num2str(theoretical_eff(ii,2),'%.3f')])
   disp(['Experimental: ',char(414),'(T1)=',num2str(experimental_eff(ii,1),'%.3f'),'   ',char(414),'(T2)=',num2str(experimental_eff(ii,2),'%.3f')])
   toc
end

%% Save results

save('steadystate_theoretical_experimental_efficiency',...
    'experimental_eff','theoretical_eff','combinations')
