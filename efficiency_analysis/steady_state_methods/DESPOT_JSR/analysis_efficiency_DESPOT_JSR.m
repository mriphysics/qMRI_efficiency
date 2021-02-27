%% Analysis of DESPOT/JSR efficiency
%David Leitao (david.leitao@kcl.ac.uk); 17-04-20

%Requires the handles to the SPGR and bSSFP steady-state signal functions 
%and its derivatives -> needs path to library folder

clearvars; close all; clc;

%% Load signal function and derivatives 
% vector u is the acquisition settings and vector p is the tissue parameters
%   u = [FA, TR, TE/RF] -> TE for SPGR sequences, RF for bSSFP sequences
%   p = [T1, T2, M0, P0, B0]

% load handles to calculate steady-state signals
load('hnd_SignalFunctions_SteadyState.mat')
signal_SPGR  =@(u,p) func_signal_SPGR (p(5),u(1),p(3),p(4),p(1),p(2),u(3),u(2));
signal_bSSFP =@(u,p) func_signal_bSSFP(p(5),u(1),p(2),p(4),u(3),p(1),p(2),u(2));

% load handles to calculate signal derivatives
load('hnd_SignalDerivatives_SteadyState.mat')
dmdT1_SPGR  =@(u,p) func_dmdT1_SPGR (p(5),u(1),p(3),p(4),p(1),p(2),u(3),u(2));
dmdT1_bSSFP =@(u,p) func_dmdT1_bSSFP(p(5),u(1),p(3),p(4),u(3),p(1),p(2),u(2));
dmdT2_SPGR  =@(u,p) func_dmdT2_SPGR (p(5),u(1),p(3),p(4),p(1),p(2),u(3),u(2));
dmdT2_bSSFP =@(u,p) func_dmdT2_bSSFP(p(5),u(1),p(3),p(4),u(3),p(1),p(2),u(2));
dmdM0_SPGR  =@(u,p) func_dmdM0_SPGR (p(5),u(1),p(4),p(1),p(2),u(3),u(2));
dmdM0_bSSFP =@(u,p) func_dmdM0_bSSFP(p(5),u(1),p(4),u(3),p(1),p(2),u(2));
dmdP0_SPGR  =@(u,p) func_dmdP0_SPGR (p(5),u(1),p(3),p(4),p(1),p(2),u(3),u(2));
dmdP0_bSSFP =@(u,p) func_dmdP0_bSSFP(p(5),u(1),p(3),p(4),u(3),p(1),p(2),u(2));
dmdB0_SPGR  =@(u,p) func_dmdB0_SPGR (p(5),u(1),p(3),p(4),p(1),p(2),u(3),u(2));
dmdB0_bSSFP =@(u,p) func_dmdB0_bSSFP(p(5),u(1),p(3),p(4),u(3),p(1),p(2),u(2));

%% Load optimised acquisition settings for DESPOT/JSR

load('opt_param_DESPOT_JSR.mat')

%% Wider range of tissue parameters for efficiency analysis

T1_list = 600:40:1200;    %Spin-lattice relaxation constant [ms]
T2_list = 40:4:100;       %Spin-spin relaxation constant [ms]
M0_list = 1;              %Equilibrium magnetisation [a.u.]
P0_list = 0;              %Constant (receiver) phase [rad]
B0_list = -100:5:100;     %Off-resonance [Hz]

nparam = numel(T1_list)*numel(T2_list)*numel(M0_list)*numel(P0_list)*numel(B0_list);
P = zeros(nparam, 5); 
cc = 1;
for i5=1:length(B0_list)
    for i4=1:numel(P0_list)
        for i3=1:numel(M0_list)
            for i2=1:length(T2_list)
                for i1=1:length(T1_list)
                    P(cc, :) = [T1_list(i1), T2_list(i2), M0_list(i3), P0_list(i4), B0_list(i5)];
                    cc = cc + 1;
                end
            end
        end
    end
end

%% Find best combination and use it for further analysis

ncomb = numel(opt_param_DESPOT_JSR.nSPGR);

[~,idx_best] = min(opt_param_DESPOT_JSR.costFunc);

uopt = [opt_param_DESPOT_JSR.acqSet{idx_best}.FAopt, ...
        opt_param_DESPOT_JSR.acqSet{idx_best}.TRopt, ...
        opt_param_DESPOT_JSR.acqSet{idx_best}.TEopt, ...
        opt_param_DESPOT_JSR.acqSet{idx_best}.RFopt];
     
% calculate T1/T2 efficiency for the wider range of tissue parameters
all_eff = zeros(nparam,2);
for nn=1:nparam
    [all_eff(nn,:), ~] = analysis_acq_set_DESPOT_JSR(uopt, P(nn,:), opt_param_DESPOT_JSR.nbSSFP(idx_best), ...
        dmdT1_SPGR,  dmdT2_SPGR,  dmdM0_SPGR,  dmdP0_SPGR,  dmdB0_SPGR, ...
        dmdT1_bSSFP, dmdT2_bSSFP, dmdM0_bSSFP, dmdP0_bSSFP, dmdB0_bSSFP,...
        signal_SPGR, signal_bSSFP);
end
best_DESPOT_JSR_T1 = reshape(all_eff(:,1), [numel(T1_list) numel(T2_list) numel(B0_list)]);
best_DESPOT_JSR_T2 = reshape(all_eff(:,2), [numel(T1_list) numel(T2_list) numel(B0_list)]);

%% Save efficiency results for best combination

save('DESPOT_JSR_eff','best_DESPOT_JSR_T1','best_DESPOT_JSR_T2')
