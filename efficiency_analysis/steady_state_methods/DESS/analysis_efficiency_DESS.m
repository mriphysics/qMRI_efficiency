%% Analysis of DESS efficiency
%David Leitao (david.leitao@kcl.ac.uk); 17-04-20

%Requires .mex file to simulate DESS steady-state signals and its 
%derivatives -> needs path to library folder

clearvars; close all; clc;

%% Load optimised acquisition settings for DESS

load('opt_param_DESS.mat')

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

%% Find best combination and use it for further analysis

ncomb = numel(opt_param_DESS.nDESS);

[~,idx_best] = min(opt_param_DESS.costFunc);

uopt = [opt_param_DESS.acqSet{idx_best}.FAopt, ...
        opt_param_DESS.acqSet{idx_best}.TRopt, ...
        opt_param_DESS.acqSet{idx_best}.TEopt];

% calculate T1/T2 efficiency for the wider range of tissue parameters
all_eff = zeros(nparam,2);
for nn=1:nparam
    [all_eff(nn,:), ~] = analysis_acq_set_DESS(uopt, P(nn,:));
end
best_DESS_T1 = reshape(all_eff(:,1), [numel(T1_list) numel(T2_list)]);
best_DESS_T2 = reshape(all_eff(:,2), [numel(T1_list) numel(T2_list)]);

%% Save efficiency results for best combination

save('DESS_eff','best_DESS_T1','best_DESS_T2')
