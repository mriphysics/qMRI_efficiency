%% Analysis of DE balanced MRF efficiency
%David Leitao (david.leitao@kcl.ac.uk); 18-04-20

clearvars; close all; clc;

%% Load all optimised acquisition settings for DE balanced MRF and merge

load('opt_param_DE_balancedMRF_ms.mat')
load('opt_param_DE_balancedMRF_ss.mat')

% merge all optimisation results
opt_param_DE_balancedMRF.N        = cat(1,opt_param_DE_balancedMRF_ms.N,opt_param_DE_balancedMRF_ss.N);
opt_param_DE_balancedMRF.costFunc = cat(1,opt_param_DE_balancedMRF_ms.costFunc,opt_param_DE_balancedMRF_ss.costFunc);
opt_param_DE_balancedMRF.acqSet   = cat(1,opt_param_DE_balancedMRF_ms.acqSet,opt_param_DE_balancedMRF_ss.acqSet);


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

%% Find best fingerprint length and use it for further analysis

% best fingerprint with N>=400 to ensure compatibility with spatial encoding
[~,idx_best] = min(opt_param_DE_balancedMRF.costFunc + 1e10*(opt_param_DE_balancedMRF.N<400));

uopt = alpha2theta(opt_param_DE_balancedMRF.acqSet{idx_best}.FAopt);
TRssfp = 5;
RF = RF_phase_cycle(opt_param_DE_balancedMRF.N(idx_best),'balanced'); %[rad]

% calculate T1/T2 efficiency for the wider range of tissue parameters
all_eff = zeros(nparam,2);
for nn=1:nparam
    [all_eff(nn,:), ~] = cppAnalysis_acq_set_DE_balancedMRF(opt_param_DE_balancedMRF.N(idx_best), ...
        uopt, RF, TRssfp, P(nn,1), P(nn,2), -P(nn,3), P(nn,5), P(nn,4));
end
best_DE_balancedMRF_T1 = reshape(all_eff(:,1), [numel(T1_list) numel(T2_list) numel(B0_list)]);
best_DE_balancedMRF_T2 = reshape(all_eff(:,2), [numel(T1_list) numel(T2_list) numel(B0_list)]);

%% Save maps of the efficiency for each method, for the best combinations

save('DE_balancedMRF_eff','best_DE_balancedMRF_T1','best_DE_balancedMRF_T2')
