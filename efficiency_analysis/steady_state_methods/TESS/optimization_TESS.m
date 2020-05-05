%% Optimisation of TESS acquisition settings 
%David Leitao (david.leitao@kcl.ac.uk); 17-04-20

%Requires .mex file to simulate TESS steady-state signals and its 
%derivatives -> needs path to library folder

clearvars; close all; clc;

%% Set tissue parameters and cost function for optimisation
% vector u is the acquisition settings and vector p is the tissue parameters
%   u = [FA, TR, TE0, dTE1] -> dTE is the time gap between consecutive echoes
%   p = [T1, T2, M0]

T1 = 781;            %Spin-lattice relaxation constant [ms]
T2 = 65;             %Spin-spin relaxation constant [ms]
M0 = 1;              %Equilibrium magnetisation [a.u.]

% set of parameters P used for optimisation
P = [T1 T2 M0];
% set cost function
func =@(u) cost_function_TESS(u, P);

%% Sequence constraints and optimisation options

TRmin = 15; %TRmin is 15MS to allow measurement of 3 echoes
TEsep = 2;  %Minimum time space between echoes and RF pulses

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
% comment lines 43-47 and replace parfor loop in line 89 by a for loop

if isempty(gcp('nocreate'))
    c = parcluster('local');
    c.NumWorkers = 4;
    parpool(c, c.NumWorkers); 
end

% Define combinations of steady-state signals to optimise
min_nTESS = 1; %minimum number of TESS steady-states (to fit 3 parameters)
max_nTESS = 6; %maximum number of TESS steady-states
comb      = min_nTESS:max_nTESS;
ncomb     = numel(comb);

% structure to save optimisation results
opt_param_TESS.nTESS    = zeros(ncomb, 1);
opt_param_TESS.costFunc = zeros(ncomb, 1);
opt_param_TESS.acqSet   = cell(ncomb, 1);

nMS = 100; %number of multi-start trials per combination
rng('default')

% optimise each combination
for nn=1:numel(comb)
    
    nTESS  = comb(nn);
    
    opt_param_TESS.nTESS(nn)  = nTESS;

    % constraints for current combination
    lb  = [zeros(nTESS,1);     TRmin*ones(nTESS,1); TEsep*ones(3*nTESS,1)];
    ub  = [pi/2*ones(nTESS,1); Inf*ones(nTESS,1);   Inf*ones(3*nTESS,1)];
    % add linear constraint TE(last) < (TR - TEsep)
    A  = zeros(1,nTESS);    b = zeros(nTESS,1);
    for ss=1:nTESS
        A(ss,nTESS+ss)   = -1; %TR(ss)
        A(ss,2*nTESS+ss) =  1; %TE0(ss)
        A(ss,3*nTESS+ss) =  1; %dTE1(ss)
        A(ss,4*nTESS+ss) =  1; %dTE2(ss)
        b(ss) = -TEsep;
    end
    % bounds used to create random initialisations (Inf cannot be used)
    lb0  = [zeros(nTESS,1);     TRmin*ones(nTESS,1); TEsep*ones(3*nTESS,1)];
    ub0  = [pi/2*ones(nTESS,1); 75*ones(nTESS,1);    50*ones(3*nTESS,1)];

    all_fval = Inf*ones(nMS, 1);
    all_uopt = zeros(nMS, length(ub));
    
    tic    
    parfor ii=1:nMS
        % random initialisation within the previously defined bounds      
        u0 = 0.8*rand(5*nTESS,1).*(ub0-lb0) + 1.1*lb0;
        
        % try random strating point; if it fails, try next one
        try
            % avoid warnings of ill-conditioning of Fisher matrix at early
            % iterations
            warning('off','all')
            [all_uopt(ii,:), all_fval(ii)] = fmincon(func, u0, A, b, [], [], lb, ub, [], options);
            warning('on','all')
        catch
            continue
        end
    end
    % extract best solution
    all_fval(all_fval<0) = Inf;
    idx_best = find(min(all_fval)==all_fval,1,'first');
    uopt = all_uopt(idx_best,:);
    
    opt_param_TESS.costFunc(nn) = all_fval(idx_best);
    
    % save optimal acquisition settings
    opt_param_TESS.acqSet{nn}.FAopt = uopt(1:nTESS);
    opt_param_TESS.acqSet{nn}.TRopt = uopt(nTESS+1:2*nTESS);
    opt_param_TESS.acqSet{nn}.TEopt = uopt(2*nTESS+1:5*nTESS);

    % display optimisation conclusion
    fprintf(1,'\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
    fprintf(1,'\nOptimisation with #TESS = %d finished.\n',nTESS)
    toc
    fprintf(1,'+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')     
end


%% Save optimisation results

% save('opt_param_TESS','opt_param_TESS')

