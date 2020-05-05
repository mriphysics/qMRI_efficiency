%% Optimisation of PLANET acquisition settings 
%David Leitao (david.leitao@kcl.ac.uk); 17-04-20

%Requires the handle to the bSSFP steady-state signal function and its 
%derivatives -> needs path to library folder

clearvars; close all; clc;

%% Load signal function and derivatives 
% vector u is the acquisition settings and vector p is the tissue parameters
%   u = [FA, TR, RF]
%   p = [T1, T2, M0, P0, B0]

% load handles to calculate steady-state signals
load('hnd_SignalFunctions_SteadyState.mat')
signal_bSSFP =@(u,p) func_signal_bSSFP(p(5),u(1),p(3),p(4),u(3),p(1),p(2),u(2));

% load handles to calculate signal derivatives
load('hnd_SignalDerivatives_SteadyState.mat')
dmdT1_bSSFP =@(u,p) func_dmdT1_bSSFP(p(5),u(1),p(3),p(4),u(3),p(1),p(2),u(2));
dmdT2_bSSFP =@(u,p) func_dmdT2_bSSFP(p(5),u(1),p(3),p(4),u(3),p(1),p(2),u(2));
dmdM0_bSSFP =@(u,p) func_dmdM0_bSSFP(p(5),u(1),p(4),u(3),p(1),p(2),u(2));
dmdP0_bSSFP =@(u,p) func_dmdP0_bSSFP(p(5),u(1),p(3),p(4),u(3),p(1),p(2),u(2));
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
% comment lines 68-72 and replace parfor loop in line 110 by a for loop

if isempty(gcp('nocreate'))
    c = parcluster('local');
    c.NumWorkers = 4;
    parpool(c, c.NumWorkers); 
end

% define combinations of bSSFP steady-state signals to optimise
min_nbSSFP = 3;  %minimum number of bSSFP steady-states (to fit 5 parameters)
max_nbSSFP = 20; %maximum number of bSSFP steady-states
comb       = min_nbSSFP:max_nbSSFP;
ncomb      = numel(comb);

% structure to save optimisation results
opt_param_PLANET.nbSSFP   = zeros(ncomb, 1);
opt_param_PLANET.costFunc = zeros(ncomb, 1);
opt_param_PLANET.acqSet   = cell(size(comb,1), 1);

nMS = 100; %number of multi-start trials per combination
rng('default')

% set cost function
func =@(u,p) cost_function_PLANET(u, p, dmdT1_bSSFP, dmdT2_bSSFP, dmdM0_bSSFP, dmdP0_bSSFP, dmdB0_bSSFP);
CostFunc = costFunc_sampling(func,P);

% optimise each combination
for nn=1:numel(comb)
    
    nbSSFP = comb(nn); 
    
    opt_param_PLANET.nbSSFP(nn) = nbSSFP;

    % constraints for current combination
    lb = [0;    TRmin;   0*ones(nbSSFP,1)];
    ub = [pi/2; Inf;     2*pi*ones(nbSSFP,1)];
    % bounds used to create random initialisations (Inf cannot be used)
    lb0 = [0;    TRmin;  0*ones(nbSSFP,1)];
    ub0 = [pi/2; 50;     2*pi*ones(nbSSFP,1)];
    
    all_fval = Inf*ones(nMS, 1);
    all_uopt = zeros(nMS, length(ub));
    
    tic    
    parfor ii=1:nMS
        % random initialisation within the previously defined bounds  
        u0 = 0.8*rand(nbSSFP+2,1).*(ub0-lb0) + 1.1*lb0;
        
        % try random strating point; if it fails, try next one
        try
            % avoid warnings of ill-conditioning of Fisher matrix at early
            % iterations
            warning('off','all')
            [all_uopt(ii,:), all_fval(ii)] = fmincon(CostFunc, u0, [], [], [], [], lb, ub, [], options);
            warning('on','all')
        catch
            continue
        end
    end
    % extract best solution
    all_fval(all_fval<0) = Inf; 
    idx_best = find(min(all_fval)==all_fval,1,'first');
    uopt = all_uopt(idx_best,:);
    
    opt_param_PLANET.costFunc(nn) = all_fval(idx_best);
    
    % save optimal acquisition settings
    opt_param_PLANET.acqSet{nn}.FAopt = uopt(1);
    opt_param_PLANET.acqSet{nn}.TRopt = uopt(2);
    opt_param_PLANET.acqSet{nn}.RFopt = uopt(3:2+nbSSFP);

    % display optimisation conclusion
    fprintf(1,'\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
    fprintf(1,'\nOptimisation with #bSSFP = %d finished.\n',nbSSFP)
    toc
    fprintf(1,'+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n') 
end


%% Save optimisation results

save('opt_param_PLANET','opt_param_PLANET')
