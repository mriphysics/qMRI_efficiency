%% Multi-start optimization of non-DE spoiled MRF acquisition settings
%David Leitao (david.leitao@kcl.ac.uk); 18-04-20

clearvars; close all; clc;

%% Set tissue parameters

T1 = 781;            %Spin-lattice relaxation constant [ms]
T2 = 65;             %Spin-spin relaxation constant [ms]
M0 = 1;              %Equilibrium magnetisation [a.u.]

% set of parameters P used for optimisation
P = [T1 T2 M0];

%% Sequence constraints and optimisation options

TRmin = 5;  %[ms]
TEmin = 2;  %[ms]
FAmax = 90; %[deg]

options = optimoptions('fmincon',...
    'SpecifyConstraintGradient',false,...
    'Algorithm','sqp',...
    'Display','none',...
    'OptimalityTolerance',1e-4,...
    'StepTolerance',1e-4,...
    'MaxIterations',1e4,...
    'MaxFunctionEvaluations',Inf,...
    'Useparallel',false);

%% Multi-start optimisation for several fingerprint lengths

% /!\ This section makes use of parallel acceleration: to disable it 
% comment lines 36-40 and replace parfor loop in line 131 by a for loop

% if isempty(gcp('nocreate'))
%     c = parcluster('local');
%     c.NumWorkers = 32;
%     parpool(c, c.NumWorkers);
% end

% define fingerprint lengths to optimise
N = [5 10 20 50 100 200]; 

% structure to save optimisation results
opt_param_nonDE_spoiledMRF_ms.N        = zeros(length(N), 1);
opt_param_nonDE_spoiledMRF_ms.costFunc = zeros(length(N), 1);
opt_param_nonDE_spoiledMRF_ms.acqSet   = cell(length(N), 1);

nMS = 100; %number of multi-start trials per fingerprint length
rng('default')

% optimise each fingerprint length
for nn=1:length(N)
     
    opt_param_nonDE_spoiledMRF_ms.N(nn) = N(nn);

    % constraints for current fingerprint length
    lb = [zeros(N(nn),1);                     TRmin*ones(N(nn),1)];
    ub = [pi; FAmax*(pi/180)*ones(N(nn)-1,1); Inf*ones(N(nn),1)];
    % bounds used to create random initialisations (Inf cannot be used)
    lb0 = [zeros(N(nn),1);                      TRmin*ones(N(nn),1)];
    ub0 = [pi; FAmax*(pi/180)*ones(N(nn)-1,1);  50*ones(N(nn),1)];
    
    all_fval = Inf*ones(nMS, 1);
    all_uopt = zeros(nMS, 2*N(nn));
    
    RF0 = zeros(N(nn), 1);      %[rad]
    TE0 = TEmin*ones(N(nn), 1); %[ms]
    
    % set cost function
    func = @(u) cppEPG_GRE_efficiency(N(nn), u(1:N(nn)), RF0, u(N(nn)+1:end), TE0, P(1), P(2));

    tic
    parfor ii=1:nMS
        % random initialisation within the lower and upper bounds   
        u0 = 0.8*rand(length(lb),1 ).*(ub0(:)-lb0(:)) + 1.1*lb0(:);
        
        % try random strating point; if it fails, try next one
        try
            % avoid warnings of ill-conditioning of Fisher matrix at early
            % iterations
            warning('off','all')
            [all_uopt(ii,:),all_fval(ii)] = fmincon(func, u0, [], [], [], [], lb, ub, [], options);
            warning('on','all')
        catch
            continue;
        end 
    end
    % extract best solution
    all_fval(all_fval<0) = Inf; all_fval(all_fval==0) = Inf;
    idx_best = find(min(all_fval)==all_fval,1,'first');
    uopt = all_uopt(idx_best,:);
    
    opt_param_nonDE_spoiledMRF_ms.costFunc(nn) = all_fval(idx_best);
    
    % save optimal acquisition settings
    opt_param_nonDE_spoiledMRF_ms.acqSet{nn}.FAopt = uopt(1:N(nn));
    opt_param_nonDE_spoiledMRF_ms.acqSet{nn}.TRopt = uopt(N(nn)+1:2*N(nn));
    
    % display optimisation conclusion
    fprintf(1,'\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
    fprintf(1,'\nOptimisation of fingerprint with length N = %d finished.\n',N(nn))
    toc
    fprintf(1,'+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
end

%% Save optimisation results

save('opt_param_nonDE_spoiledMRF_ms','opt_param_nonDE_spoiledMRF_ms')

