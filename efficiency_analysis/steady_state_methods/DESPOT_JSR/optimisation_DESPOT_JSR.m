%% Optimisation of DESPOT/JSR acquisition settings 
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
signal_bSSFP =@(u,p) func_signal_bSSFP(p(5),u(1),p(3),p(4),u(3),p(1),p(2),u(2));

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

% set of parameters P used for optimisation
T1_list = 781;            %Spin-lattice relaxation constant [ms]
T2_list = 65;             %Spin-spin relaxation constant [ms]
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

%% Sequence constraints and optimisation options

TRmin = 5;
TEmin = 2; %for SPGR only
TEsep = 2; %for SPGR only

options = optimoptions('fmincon',...
    'SpecifyConstraintGradient',false,...
    'Algorithm','sqp',...
    'MaxFunctionEvaluations',Inf,...
    'Display','none',...
    'OptimalityTolerance',1e-4,...
    'StepTolerance',1e-4,...
    'MaxIterations',1e4,...
    'Useparallel',false);

%% Optimisation for several combinations of steady-state sequences

% /!\ This section makes use of parallel acceleration: to disable it 
% comment lines 76-80 and replace parfor loop in line 131 by a for loop

if isempty(gcp('nocreate'))
    c = parcluster('local');
    c.NumWorkers = 4;
    parpool(c, c.NumWorkers); 
end

% define combinations of steady-state signals to optimise
min_nSS = 3; %minimum number of steady-states (to fit 5 parameters)
max_nSS = 8; %maximum number of steady-states
comb    = gen_comb_steady_states(max_nSS,min_nSS); %generate all possible combinations    
comb(comb(:,2)==0,:) = []; %remove unfeasible combinations with 0 bSSFP
ncomb = size(comb,1);

% structure to save optimisation results
opt_param_DESPOT_JSR.nSPGR    = zeros(ncomb, 1);
opt_param_DESPOT_JSR.nbSSFP   = zeros(ncomb, 1);
opt_param_DESPOT_JSR.costFunc = zeros(ncomb, 1);
opt_param_DESPOT_JSR.acqSet   = cell(ncomb, 1);

nMS = 100; %number of multi-start trials per combination
rng('default')

% optimise each combination
for nn=1:ncomb
    
    nSPGR  = comb(nn,1);
    nbSSFP = comb(nn,2);
    nSS = nSPGR + nbSSFP;
    
    opt_param_DESPOT_JSR.nSPGR(nn)  = nSPGR;
    opt_param_DESPOT_JSR.nbSSFP(nn) = nbSSFP;

    % constraints for current combination
    lb = [0*(pi/180)*ones(nSS,1); TRmin*ones(nSS,1); TEmin*ones(nSPGR,1);   0*ones(nbSSFP,1)];
    ub = [pi/2*ones(nSS,1);       Inf*ones(nSS,1);   Inf*ones(nSPGR,1);     2*pi*ones(nbSSFP,1)];
    A  = zeros(nSPGR,3*nSS);	b = zeros(nSPGR,1);
    for ss=1:nSPGR
        A(ss,2*nSS+ss) = 1; %TE(ss)
        A(ss,nSS+ss) = -1;  %TR(ss)
        b(ss) = -TEsep;
    end
    % bounds used to create random initialisations (Inf cannot be used)
    lb0 = [0*(pi/180)*ones(nSS,1); TRmin*ones(nSS,1); TEmin*ones(nSPGR,1);      0*ones(nbSSFP,1)];
    ub0 = [pi/2*ones(nSS,1);       50*ones(nSS,1);    (50-TEsep)*ones(nSPGR,1); 2*pi*ones(nbSSFP,1)];

    all_fval = Inf*ones(nMS, 1);
    all_uopt = zeros(nMS, length(ub));
    
    % set cost function
    func =@(u,p) cost_function_DESPOT_JSR(u, p, nbSSFP, ...
        dmdT1_SPGR,  dmdT2_SPGR,  dmdM0_SPGR,  dmdP0_SPGR,  dmdB0_SPGR, ...
        dmdT1_bSSFP, dmdT2_bSSFP, dmdM0_bSSFP, dmdP0_bSSFP, dmdB0_bSSFP);
    CostFunc = costFunc_sampling(func,P);
    
    tic    
    parfor ii=1:nMS
        % random initialisation within the previously defined bounds       
        u0 = 0.8*rand(3*nSS,1).*(ub0-lb0) + 1.1*lb0;
        
        % try random strating point; if it fails, try next one
        try
            % avoid warnings of ill-conditioning of Fisher matrix at early
            % iterations
            warning('off','all')
            [all_uopt(ii,:), all_fval(ii)] = fmincon(CostFunc, u0, A, b, [], [], lb, ub, [], options);
            warning('on','all')
        catch
            continue
        end
    end
    % extract best solution
    all_fval(all_fval<0) = Inf;
    idx_best = find(min(all_fval)==all_fval,1,'first');
    uopt = all_uopt(idx_best,:);
    
    opt_param_DESPOT_JSR.costFunc(nn) = all_fval(idx_best);
    
    % save optimal acquisition settings
    opt_param_DESPOT_JSR.acqSet{nn}.FAopt = uopt(1:nSS);
    opt_param_DESPOT_JSR.acqSet{nn}.TRopt = uopt(nSS+1:2*nSS); 
    opt_param_DESPOT_JSR.acqSet{nn}.TEopt = uopt(2*nSS+1:2*nSS+nSPGR);
    opt_param_DESPOT_JSR.acqSet{nn}.RFopt = uopt(2*nSS+nSPGR+1:3*nSS);
    
    % display optimisation conclusion
    fprintf(1,'\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
    fprintf(1,'\nOptimisation with {#SPGR:#bSSFP} = {%d:%d} finished.\n',nSPGR,nbSSFP)
    toc
    fprintf(1,'+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
end


%% Save optimisation results

save('opt_param_DESPOT_JSR','opt_param_DESPOT_JSR')
