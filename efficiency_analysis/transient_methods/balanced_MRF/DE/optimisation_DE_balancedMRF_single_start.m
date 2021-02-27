%% Single-start optimisation of DE balanced MRF acquisition settings
%David Leitao (david.leitao@kcl.ac.uk); 18-04-20

clearvars; close all; clc;

%% Set tissue parameters and initial acquisition settings

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

% initial acquisition settings; adapted from Ma et al. 2013 (DOI:10.1038/nature11971)
theta0 = deg2rad([(10+50*sin(2*pi*(1:250)/500)), (10+20*sin(2*pi*(1:250)/500)), ...
    (10+50*sin(2*pi*(1:250)/500)), (10+20*sin(2*pi*(1:250)/500))])';

%% Sequence constraints and optimisation options

TRssfp = 5;  %[ms]
FAmax  = 90; %[deg]

options = optimoptions('fmincon',...
    'SpecifyConstraintGradient',false,...
    'Algorithm','sqp',...
    'Display','none',...
    'OptimalityTolerance',1e-4,...
    'StepTolerance',1e-4,...
    'MaxIterations',1e4,...
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
N = [500 600 700 800 900 1000];

% structure to save optimisation results
opt_param_DE_balancedMRF_ss.N        = zeros(length(N), 1);
opt_param_DE_balancedMRF_ss.costFunc = zeros(length(N), 1);
opt_param_DE_balancedMRF_ss.acqSet   = cell(length(N), 1);

rng('default')

% optimise each fingerprint length
for nn=1:length(N)
     
    opt_param_DE_balancedMRF_ss.N(nn) = N(nn);
       
    % constraints for current fingerprint length
    lb = zeros(N(nn),1);
    ub = [0; (FAmax/2)*(pi/180)*ones(N(nn)-1,1)];
    % first derivative constraint
    a1 = 5; 
    A1 = zeros(N(nn),N(nn));
    A1(1+N(nn)+(0:N(nn)-2)*(N(nn)+1)) = 1; A1(N(nn)) = 1;
    A1(1+(0:N(nn)-1)*(N(nn)+1)) = -1;
    b1 = deg2rad(a1)*ones(2*N(nn),1);
    % second derivative constraint
    a2 = 0.5;
    A2 = zeros(N(nn),N(nn));
    A2(1:N(nn)+1:N(nn)*N(nn)) = 1;
    A2(N(nn)+1:N(nn)+1:N(nn)*N(nn)) = -2; A2(N(nn)) = -2;
    A2(2*N(nn)+1:N(nn)+1:N(nn)*N(nn)) = 1; A2(N(nn)-1) = 1; A2(2*N(nn)) = 1;
    b2 = deg2rad(a2)*ones(2*N(nn),1);

    A = [A1; -A1; A2; -A2];
    b = [b1; b2];
    
    RF = RF_phase_cycle(N(nn),'balanced'); %[rad]    
    
    % set cost function
    func =@(u,p) cppDE_balancedMRF_efficiency(N(nn), u, RF, TRssfp, ...
        p(1), p(2), -p(3), p(5), p(4)); % start with magnetisation inverted: -p(3)
    CostFunc = costFunc_sampling(func,P);

    u0 = theta0(1:N(nn));
    
    tic
    [uopt, fval] = fmincon(CostFunc, u0, A, b, [], [], lb, ub, [], options);

    opt_param_DE_balancedMRF_ss.costFunc(nn) = fval;
    
    % save optimal acquisition settings
    opt_param_DE_balancedMRF_ss.acqSet{nn}.FAopt = theta2alpha(uopt);
    
    % display optimisation conclusion
    fprintf(1,'\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
    fprintf(1,'\nOptimisation of fingerprint with length N = %d finished.\n',N(nn))
    toc
    fprintf(1,'+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
end

%% Save optimisation results

save('opt_param_DE_balancedMRF_ss','opt_param_DE_balancedMRF_ss')

