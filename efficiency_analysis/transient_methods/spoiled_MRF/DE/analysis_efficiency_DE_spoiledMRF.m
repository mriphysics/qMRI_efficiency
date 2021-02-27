%% Analysis of DE spoiled MRF efficiency
%David Leitao (david.leitao@kcl.ac.uk); 18-04-20

clearvars; close all; clc;

%% Load all optimised acquisition settings for DE spoiled MRF and merge

load('opt_param_DE_spoiledMRF_ms.mat')
load('opt_param_DE_spoiledMRF_ss.mat')

% merge all optimisation results
opt_param_DE_spoiledMRF.N        = cat(1,opt_param_DE_spoiledMRF_ms.N,opt_param_DE_spoiledMRF_ss.N);
opt_param_DE_spoiledMRF.costFunc = cat(1,opt_param_DE_spoiledMRF_ms.costFunc,opt_param_DE_spoiledMRF_ss.costFunc);
opt_param_DE_spoiledMRF.acqSet   = cat(1,opt_param_DE_spoiledMRF_ms.acqSet,opt_param_DE_spoiledMRF_ss.acqSet);


%% Wider range of tissue parameters for efficiency analysis

T1_list = 600:40:1200;    %Spin-lattice relaxation constant [ms]
T2_list = 40:4:100;       %Spin-spin relaxation constant [ms]
M0_list = 1;              %Equilibrium magnetisation [a.u.]

nparam = numel(T1_list)*numel(T2_list)*numel(M0_list);
P = zeros(nparam, 3); 
cc = 1;
for i3=1:numel(M0_list)
    for i2=1:length(T2_list)
        for i1=1:length(T1_list)
            P(cc, :) = [T1_list(i1), T2_list(i2), M0_list(i3)];
            cc = cc + 1;
        end
    end
end

%% Find best fingerprint length and use it for further analysis

% best fingerprint with N>=400 to ensure compatibility with spatial encoding
[~,idx_best] = min(opt_param_DE_spoiledMRF.costFunc + 1e10*(opt_param_DE_spoiledMRF.N<400));

RF0 = zeros(opt_param_DE_spoiledMRF.N(idx_best), 1);      %[rad]
TEmin = 2;                                                %[ms]
TE0 = TEmin*ones(opt_param_DE_spoiledMRF.N(idx_best), 1); %[ms]

% calculate T1/T2 efficiency for the wider range of tissue parameters
all_eff = zeros(nparam,2);
for nn=1:nparam         %cppAnalysis_acq_set_DE_spoiledMRF
    [all_eff(nn,:), ~] = cppAnalysis_acq_set_DE_spoiledMRF(opt_param_DE_spoiledMRF.N(idx_best), ...
        opt_param_DE_spoiledMRF.acqSet{idx_best}.FAopt, RF0, ...
        opt_param_DE_spoiledMRF.acqSet{idx_best}.TRopt, TE0, P(nn,1), P(nn,2));
end
best_DE_spoiledMRF_T1 = reshape(all_eff(:,1), [numel(T1_list) numel(T2_list)]);
best_DE_spoiledMRF_T2 = reshape(all_eff(:,2), [numel(T1_list) numel(T2_list)]);

%% Save maps of the efficiency for each method, for the best combinations

save('DE_spoiledMRF_eff','best_DE_spoiledMRF_T1','best_DE_spoiledMRF_T2')
