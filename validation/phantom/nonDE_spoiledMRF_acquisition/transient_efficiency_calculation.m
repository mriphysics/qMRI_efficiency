%% Calculation and comparison of experimental and theoretical efficiency
%David Leitao (david.leitao@kcl.ac.uk); 24-04-20

clearvars; close all; clc;

%% Define acquisition settings used  

Nmin = 1;
Nmax = 800;
N    = Nmax-Nmin+1;

FA = load('Jiang_FA_FISP.txt'); 
FA = deg2rad(FA(Nmin:Nmax));
RF = zeros(N,1);
TR = 7.0*ones(N,1);
TE = 2.3*ones(N,1);

%% Load data

% load gold standard parameters 
load('transient_gold_standard_parameters.mat')

% load data
load('MRF_posproc.mat')

x_pnt = 76:84; %spatial points used to boost amount of data
npos  = length(x_pnt);

% remove first dynamic and extract only data points taht will be used
MRF_data = squeeze(MRF_all(x_pnt, Nmin:Nmax, 2:end));
MRF_data = permute(MRF_data,[2 3 1]);

clear MRF_all %release memory

%% Combinations of fingerprints segments

ncomb = 0;
nseg  = 8; 
for ii=1:nseg
    ncomb = ncomb + size(combnk(1:nseg,ii),1);
end

combList = cell(ncomb,1);
w        = zeros(ncomb, N);

cc = 1;
for ii=1:nseg
    aux_comb = combnk(1:nseg,ii);
    for kk=1:size(aux_comb,1)
        combList{cc} = aux_comb(kk,:);
        cc = cc + 1; 
    end
end

for ii=1:ncomb
    for kk=1:numel(combList{ii})
        w(ii, (1+(combList{ii}(kk)-1)*(N/nseg)):(combList{ii}(kk)*(N/nseg))) = true;
    end
end
w = logical(w);

%% Calculate experimental efficiency

% /!\ This section makes use of parallel acceleration: to disable it 
% comment lines 67-71 and replace parfor loop in line 90 by a for loop

if isempty(gcp('nocreate'))
    c = parcluster('local');
    c.NumWorkers = 4;
    parpool(c, c.NumWorkers);
end

ndyn = size(MRF_data,2);
            
T1_est = zeros(ncomb, ndyn, npos);
T2_est = zeros(ncomb, ndyn, npos);
M0_est = zeros(ncomb, ndyn, npos);

all_experimental_eff = zeros(ncomb, 2, npos);
all_theoretical_eff  = zeros(ncomb, 2, npos);

for nn=1:npos

    x0 = [T1_gs(x_pnt(nn)) T2_gs(x_pnt(nn)) M0_gs(x_pnt(nn))];

    all_signals = abs(squeeze((MRF_data(:,:,nn))));

    for ii=1:ncomb
        f =@(x) signals_combination(N, B1_gs(x_pnt(nn)).*FA, RF, TR, TE, x(1), x(2), x(3), double(w(ii,:)));
        parfor kk=1:ndyn
            meas_signal = all_signals(w(ii,:), kk);

            Fmin = @(x) (f(x) - meas_signal);

            lb = [0   0   0]; 
            ub = [Inf Inf Inf];
            opts =  optimoptions('lsqnonlin',...
                'Display','off',...
                'MaxFunctionEvaluations',Inf,...
                'MaxIterations',1e4,...
                'FunctionTolerance',1e-8,...
                'StepTolerance',1e-6);
            [x_est,res_norm] = lsqnonlin(Fmin,x0,lb,ub,opts);

            T1_est(ii, kk, nn) = x_est(1);
            T2_est(ii, kk, nn) = x_est(2);
            M0_est(ii, kk, nn) = x_est(3);        
        end

        T1_std = std(T1_est(ii,:,nn));
        T2_std = std(T2_est(ii,:,nn));
       
        all_experimental_eff(ii,1,nn) = T1_gs(x_pnt(nn)) * noise_std(x_pnt(nn),channel) / ...
            (T1_std * M0_gs(x_pnt(nn)) * sqrt(sum(TR(w(ii,:)))*1e-3));
        all_experimental_eff(ii,2,nn) = T2_gs(x_pnt(nn)) * noise_std(x_pnt(nn),channel) / ...
            (T2_std * M0_gs(x_pnt(nn)) * sqrt(sum(TR(w(ii,:)))*1e-3));
        
        %cpp function returns efficiency calculated dusing all TRs, hence
        %the rescaling
        all_theoretical_eff(ii, :, nn) = ...
            cppEPG_GRE_weighted_efficiency(N, double(w(ii,:)), B1_gs(x_pnt(nn)).*FA, RF, ...
            TR, TE, T1_gs(x_pnt(nn)), T2_gs(x_pnt(nn))) * sqrt(sum(TR)*1e-3)/sqrt(sum(TR(w(ii,:)))*1e-3);

       disp(' ')
       disp([num2str(ii),'/',num2str(ncomb)])
       disp(['Theoretical:  ',char(414),'(T1)=',num2str(sum(all_theoretical_eff(ii,1,:),3)/nn,'%.3f'),'   ',char(414),'(T2)=',num2str(sum(all_theoretical_eff(ii,2,:),3)/nn,'%.3f')])
       disp(['Experimental: ',char(414),'(T1)=',num2str(sum(all_experimental_eff(ii,1,:),3)/nn,'%.3f'),'   ',char(414),'(T2)=',num2str(sum(all_experimental_eff(ii,2,:),3)/nn,'%.3f')])
    end
end

%% Save results

save('transient_theoretical_experimental_efficiency',...
     'all_experimental_eff','all_theoretical_eff','w')
