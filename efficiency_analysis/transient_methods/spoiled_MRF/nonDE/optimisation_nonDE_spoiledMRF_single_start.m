%% Single-start optimisation of non-DE spoiled MRF acquisition settings
%David Leitao (david.leitao@kcl.ac.uk); 18-04-20

clearvars; close all; clc;

%% Set tissue parameters and initial acquisition settings

T1 = 781;            %Spin-lattice relaxation constant [ms]
T2 = 65;             %Spin-spin relaxation constant [ms]
M0 = 1;              %Equilibrium magnetisation [a.u.]

% set of parameters P used for optimisation
P = [T1 T2 M0];

% initial acquisition settings; Jiang et al. 2015 (DOI:10.1002/mrm.25559)
FA0 = load('Jiang_FA_FISP.txt'); 
FA0 = deg2rad([FA0(:); FA0(:)]);
TR0 = load('Jiang_TR_FISP.txt');
TR0 = [TR0(:); TR0(:)];

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
    'MaxIterations',4,...
    'MaxFunctionEvaluations',Inf,...
    'Useparallel',true);

%% Single-start optimisation for several fingerprint lengths

% /!\ This section makes use of parallel acceleration: to disable it 
% set 'Useparallel' to false in line 48

% if isempty(gcp('nocreate'))
%     c = parcluster('local');
%     c.NumWorkers = 32;
%     parpool(c, c.NumWorkers);
% end

% define fingerprint lengths to optimise
N = [300 400 500 600 700 800 900 1000 1200 1400 1600];

% structure to save optimisation results
opt_param_nonDE_spoiledMRF_ss.N        = zeros(length(N), 1);
opt_param_nonDE_spoiledMRF_ss.costFunc = zeros(length(N), 1);
opt_param_nonDE_spoiledMRF_ss.acqSet   = cell(length(N), 1);

rng('default')

% optimise each fingerprint length
for nn=1:length(N)
     
    opt_param_nonDE_spoiledMRF_ss.N(nn) = N(nn);

    % constraints for current fingerprint length
    lb = [zeros(N(nn),1);                     TRmin*ones(N(nn),1)];
    ub = [pi; FAmax*(pi/180)*ones(N(nn)-1,1); Inf*ones(N(nn),1)];
    
    RF0 = zeros(N(nn), 1);      %[rad]
    TE0 = TEmin*ones(N(nn), 1); %[ms]
    
    % set cost function
    func = @(u) cppEPG_GRE_efficiency(N(nn), u(1:N(nn)), RF0, u(N(nn)+1:end), TE0, T1, T2);
    
    u0 = [pi; FA0(1:N(nn)-1); TR0(1:N(nn))];

    tic
    [uopt,fval] = fmincon(func, u0, [], [], [], [], lb, ub, [], options);

    opt_param_nonDE_spoiledMRF_ss.costFunc(nn) = fval;
    
    % save optimal acquisition settings
    opt_param_nonDE_spoiledMRF_ss.acqSet{nn}.FAopt = uopt(1:N(nn));
    opt_param_nonDE_spoiledMRF_ss.acqSet{nn}.TRopt = uopt(N(nn)+1:2*N(nn));
    
    % display optimisation conclusion
    fprintf(1,'\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
    fprintf(1,'\nOptimisation of fingerprint with length N = %d finished.\n',N(nn))
    toc
    fprintf(1,'+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n') 
end

%% Save optimisation results

save('opt_param_nonDE_spoiledMRF_ss','opt_param_nonDE_spoiledMRF_ss')

