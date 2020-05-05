%% Estimate gold standard tissue parameters
%David Leitao (david.leitao@kcl.ac.uk); 24-04-20

clearvars; close all; clc;

%% Load data and define acquisition settings used

N = 800;

FA = load('Jiang_FA_FISP.txt'); 
FA = deg2rad(FA(1:N));

RF = zeros(N,1);
TR = 7.0*ones(N,1);
TE = 2.3*ones(N,1);

load('MRF_posproc.mat')

%% Fit gold standard values of [T1,T2,M0,B1]

n_pts = size(MRF_all, 1);
            
T1_gs = zeros(n_pts, 1);
T2_gs = zeros(n_pts, 1);
M0_gs = zeros(n_pts, 1);
B1_gs = zeros(n_pts, 1);

for nn=1:n_pts
    %    [T1   T2   M0    B1]
    x0 = [1200 120  1600  1]; %<- starting point
    lb = [0    0    0     0]; 
    ub = [Inf  Inf  Inf   Inf];
    
    f =@(x) signals_combination(N, x(4).*FA, RF, TR, TE, x(1), x(2), x(3), ones(N,1));
    mean_signal = abs(mean(squeeze(MRF_all(nn,:,:)),2));
    Fmin = @(x) (f(x) - mean_signal);
   
    opts =  optimoptions('lsqnonlin',...
        'Display','none',...
        'MaxFunctionEvaluations',Inf,...
        'MaxIterations',1e4,...
        'FunctionTolerance',1e-8);

    [x_est,res_norm] = lsqnonlin(Fmin,x0,lb,ub,opts);

    T1_gs(nn) = x_est(1);
    T2_gs(nn) = x_est(2);
    M0_gs(nn) = x_est(3);    
    B1_gs(nn) = x_est(4); 
    
    disp([num2str(nn),'/',num2str(n_pts)])
end

%% Save gold standard values

save('transient_gold_standard_parameters','T1_gs','T2_gs','M0_gs','B1_gs','noise_std')

%% Plot fitted values

figure; 
set(gcf,'Units','normalized','Outerposition',[0 0 1 1],'Color','w')
subplot(2,2,1) 
plot(T1_gs,'Linewidth',2)
xlim([1 120]); ylim([0 2000]); grid minor;
ylabel('T_1 (ms)'); xlabel('x position (FE direction)'); title('T_1')
set(gca,'Fontsize',14)

subplot(2,2,2) 
plot(T2_gs,'Linewidth',2)
xlim([1 120]); ylim([0 200]); grid minor;
ylabel('T_2 (ms)'); xlabel('x position (FE direction)'); title('T_2')
set(gca,'Fontsize',14)

subplot(2,2,3) 
plot(M0_gs,'Linewidth',2)
xlim([1 120]); ylim([0 2000]); grid minor;
ylabel('M_0 (a.u.)'); xlabel('x position (FE direction)'); title('M_0')
set(gca,'Fontsize',14)

subplot(2,2,4) 
plot(B1_gs,'Linewidth',2)
xlim([1 120]); ylim([0.8 1.1]); grid minor;
ylabel('B_1 (a.u.)'); xlabel('x position (FE direction)'); title('B_1')
set(gca,'Fontsize',14)
