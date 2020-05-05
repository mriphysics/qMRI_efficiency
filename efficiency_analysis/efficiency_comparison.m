%% Comparison of the efficiency distribution across analysed methods
%David Leitao (david.leitao@kcl.ac.uk); 18-04-20

clearvars; close all; clc;

%% Wider range of tissue parameters for efficiency comparison

T1_list = 600:40:1200;    %Spin-lattice relaxation constant [ms]
T2_list = 40:4:100;       %Spin-spin relaxation constant [ms]
B0_list = -100:5:100;     %Off-resonance [Hz]

%% Load optimised results for each method

% steady-state methods
load('DESPOT_JSR_eff.mat')
load('PLANET_eff.mat')
load('DESS_eff.mat')
load('TESS_eff.mat')

% transient methods
load('nonDE_spoiledMRF_eff.mat')
load('DE_spoiledMRF_eff.mat') 
load('nonDE_balancedMRF_eff.mat')
load('DE_balancedMRF_eff.mat') 

%% Efficiency dispersion as a function of B0 and {T1,T2}

nT1 = numel(T1_list);
nT2 = numel(T2_list);
nB0 = numel(B0_list);

% efficiency distribution as a function of B0 and averaged over {T1,T2}
B0_DESPOT_JSR.T1 = squeeze(mean(mean(best_DESPOT_JSR_T1,1),2));
B0_DESPOT_JSR.T2 = squeeze(mean(mean(best_DESPOT_JSR_T2,1),2));
B0_PLANET.T1     = squeeze(mean(mean(best_PLANET_T1,1),2));
B0_PLANET.T2     = squeeze(mean(mean(best_PLANET_T2,1),2));
B0_DESS.T1       = squeeze(mean(mean(best_DESS_T1,1),2));
B0_DESS.T2       = squeeze(mean(mean(best_DESS_T2,1),2));
B0_TESS.T1       = squeeze(mean(mean(best_TESS_T1,1),2));
B0_TESS.T2       = squeeze(mean(mean(best_TESS_T2,1),2));
B0_DE_balancedMRF.T1    = squeeze(mean(mean(mean(best_DE_balancedMRF_T1,1),2),4));
B0_DE_balancedMRF.T2    = squeeze(mean(mean(mean(best_DE_balancedMRF_T2,1),2),4));
B0_nonDE_balancedMRF.T1 = squeeze(mean(mean(mean(best_nonDE_balancedMRF_T1,1),2),4));
B0_nonDE_balancedMRF.T2 = squeeze(mean(mean(mean(best_nonDE_balancedMRF_T2,1),2),4));
B0_DE_spoiledMRF.T1     = squeeze(mean(mean(mean(best_DE_spoiledMRF_T1,1),2),3));
B0_DE_spoiledMRF.T2     = squeeze(mean(mean(mean(best_DE_spoiledMRF_T2,1),2),3));
B0_nonDE_spoiledMRF.T1  = squeeze(mean(mean(mean(best_nonDE_spoiledMRF_T1,1),2),3));
B0_nonDE_spoiledMRF.T2  = squeeze(mean(mean(mean(best_nonDE_spoiledMRF_T2,1),2),3));

% efficiency distribution as a function of {T1,T2} and averaged over B0
T1T2_DESPOT_JSR.T1 = mean(best_DESPOT_JSR_T1,3);
T1T2_DESPOT_JSR.T2 = mean(best_DESPOT_JSR_T2,3);
T1T2_PLANET.T1     = mean(best_PLANET_T1,3);
T1T2_PLANET.T2     = mean(best_PLANET_T2,3);
T1T2_DESS.T1       = best_DESS_T1;
T1T2_DESS.T2       = best_DESS_T2;
T1T2_TESS.T1       = best_TESS_T1;
T1T2_TESS.T2       = best_TESS_T2;
T1T2_DE_balanced.T1    = mean(best_DE_balancedMRF_T1,3);
T1T2_DE_balanced.T2    = mean(best_DE_balancedMRF_T2,3);
T1T2_nonDE_balanced.T1 = mean(best_nonDE_balancedMRF_T1,3);
T1T2_nonDE_balanced.T2 = mean(best_nonDE_balancedMRF_T2,3);
T1T2_DE_spoiled.T1     = best_DE_spoiledMRF_T1;
T1T2_DE_spoiled.T2     = best_DE_spoiledMRF_T2;
T1T2_nonDE_spoiled.T1  = best_nonDE_spoiledMRF_T1;
T1T2_nonDE_spoiled.T2  = best_nonDE_spoiledMRF_T2;

%% Plot efficiency distribution for all analysed qMRI methods

FntSz = 14;
c = lines(2);
c_sort = [c(1,:); c(1,:); c(1,:); c(1,:); c(2,:); c(2,:); c(2,:); c(2,:)];

figure;
set(gcf,'Units','normalized','Outerposition',[0.05 0 0.9 1],'Color','w');


subplot(2,2,1)
boxplot(cat(1, T1T2_DESPOT_JSR.T1(:), T1T2_PLANET.T1(:), ...
    T1T2_DESS.T1(:), T1T2_TESS.T1(:), T1T2_nonDE_spoiled.T1(:), ...
    T1T2_DE_spoiled.T1(:), T1T2_nonDE_balanced.T1(:), T1T2_DE_balanced.T1(:)),...
    cat(1,repmat('     DESPOT/JSR',[nT1*nT2 1]), ...
          repmat('         PLANET',[nT1*nT2 1]), ...
          repmat('           DESS',[nT1*nT2 1]), ...
          repmat('           TESS',[nT1*nT2 1]), ...
          repmat(' non-DE spoiled',[nT1*nT2 1]), ...
          repmat('     DE spoiled',[nT1*nT2 1]), ...
          repmat('non-DE balanced',[nT1*nT2 1]), ...
          repmat('    DE balanced',[nT1*nT2 1])), ...
    'Colors',c_sort)
set(gca,'XTickLabelRotation',20); 
set(findobj(gca,'type','line'),'linew',1.5)
ylim([0 1.3]); yticks(0:0.3:1.2)
ylabel('\eta'); set(get(gca,'ylabel'),'rotation',0); 
set(gca,'Fontsize',FntSz)
set(get(gca,'ylabel'),'units','normalized','position',[-0.19 0.53],'Fontsize',FntSz+6)
text(-0.175,0.555,'(T_1)','Fontsize',FntSz+2,'units','normalized')
text(-0.2,0.4,'(s^{-1/2})','Fontsize',FntSz,'units','normalized')
text(-0.25,1.05,'a','Fontsize',FntSz+8,'units','normalized','Fontweight','bold')
text(1.66,1.1,'                ','Fontsize',FntSz,...
    'HorizontalAlignment','center','Edgecolor',[0 0 0])
text(1.33,1.11,'\color{red}+','Fontsize',FntSz,...
    'HorizontalAlignment','right')
text(1.3,1.1,'\color{black} Outlier','Fontsize',FntSz,...
    'HorizontalAlignment','left')


subplot(2,2,2)
boxplot(cat(1, T1T2_DESPOT_JSR.T2(:), T1T2_PLANET.T2(:), ...
    T1T2_DESS.T2(:), T1T2_TESS.T2(:), T1T2_nonDE_spoiled.T2(:), ...
    T1T2_DE_spoiled.T2(:), T1T2_nonDE_balanced.T2(:), T1T2_DE_balanced.T2(:)),...
    cat(1,repmat('     DESPOT/JSR',[nT1*nT2 1]), ...
          repmat('         PLANET',[nT1*nT2 1]), ...
          repmat('           DESS',[nT1*nT2 1]), ...
          repmat('           TESS',[nT1*nT2 1]), ...
          repmat(' non-DE spoiled',[nT1*nT2 1]), ...
          repmat('     DE spoiled',[nT1*nT2 1]), ...
          repmat('non-DE balanced',[nT1*nT2 1]), ...
          repmat('    DE balanced',[nT1*nT2 1])), ...
    'Colors',c_sort)
set(gca,'XTickLabelRotation',20);
set(findobj(gca,'type','line'),'linew',1.5)
ylim([0 1.3]); yticks(0:0.3:1.2)
ylabel('\eta'); set(get(gca,'ylabel'),'rotation',0); 
set(gca,'Fontsize',FntSz)
set(get(gca,'ylabel'),'units','normalized','position',[-0.19 0.53],'Fontsize',FntSz+6)
text(-0.175,0.555,'(T_2)','Fontsize',FntSz+2,'units','normalized')
text(-0.2,0.4,'(s^{-1/2})','Fontsize',FntSz,'units','normalized')
text(-0.25,1.05,'b','Fontsize',FntSz+8,'units','normalized','Fontweight','bold')
text(1.66,1.1,'                ','Fontsize',FntSz,...
    'HorizontalAlignment','center','Edgecolor',[0 0 0])
text(1.33,1.11,'\color{red}+','Fontsize',FntSz,...
    'HorizontalAlignment','right')
text(1.3,1.1,'\color{black} Outlier','Fontsize',FntSz,...
    'HorizontalAlignment','left')


subplot(2,2,3)
boxplot(cat(1, B0_DESPOT_JSR.T1(:), B0_PLANET.T1(:), ...
    B0_DESS.T1(:), B0_TESS.T1(:), B0_nonDE_spoiledMRF.T1(:), ...
    B0_DE_spoiledMRF.T1(:), B0_nonDE_balancedMRF.T1(:), B0_DE_balancedMRF.T1(:)),...
    cat(1,repmat('     DESPOT/JSR',[nB0 1]), ...
          repmat('         PLANET',[nB0 1]), ...
          repmat('           DESS',[1   1]), ...
          repmat('           TESS',[1   1]), ...
          repmat(' non-DE spoiled',[1   1]), ...
          repmat('     DE spoiled',[1   1]), ...
          repmat('non-DE balanced',[nB0 1]), ...
          repmat('    DE balanced',[nB0 1])), ...
    'Colors',c_sort)
set(gca,'XTickLabelRotation',20); 
set(findobj(gca,'type','line'),'linew',1.5)
ylim([0 1.3]); yticks(0:0.3:1.2)
ylabel('\eta'); set(get(gca,'ylabel'),'rotation',0); 
set(gca,'Fontsize',FntSz)
set(get(gca,'ylabel'),'units','normalized','position',[-0.19 0.53],'Fontsize',FntSz+6)
text(-0.175,0.555,'(T_1)','Fontsize',FntSz+2,'units','normalized')
text(-0.2,0.4,'(s^{-1/2})','Fontsize',FntSz,'units','normalized')
text(-0.25,1.05,'c','Fontsize',FntSz+8,'units','normalized','Fontweight','bold')
text(1.66,1.1,'                ','Fontsize',FntSz,...
    'HorizontalAlignment','center','Edgecolor',[0 0 0])
text(1.33,1.11,'\color{red}+','Fontsize',FntSz,...
    'HorizontalAlignment','right')
text(1.3,1.1,'\color{black} Outlier','Fontsize',FntSz,...
    'HorizontalAlignment','left')


subplot(2,2,4)
boxplot(cat(1, B0_DESPOT_JSR.T2(:), B0_PLANET.T2(:), ...
    B0_DESS.T2(:), B0_TESS.T2(:), B0_nonDE_spoiledMRF.T2(:), ...
    B0_DE_spoiledMRF.T2(:), B0_nonDE_balancedMRF.T2(:), B0_DE_balancedMRF.T2(:)),...
    cat(1,repmat('     DESPOT/JSR',[nB0 1]), ...
          repmat('         PLANET',[nB0 1]), ...
          repmat('           DESS',[1   1]), ...
          repmat('           TESS',[1   1]), ...
          repmat(' non-DE spoiled',[1   1]), ...
          repmat('     DE spoiled',[1   1]), ...
          repmat('non-DE balanced',[nB0 1]), ...
          repmat('    DE balanced',[nB0 1])), ...
    'Colors',c_sort)
set(gca,'XTickLabelRotation',20);
set(findobj(gca,'type','line'),'linew',1.5)
ylim([0 1.3]); yticks(0:0.3:1.2)
ylabel('\eta'); set(get(gca,'ylabel'),'rotation',0); 
set(gca,'Fontsize',FntSz)
set(get(gca,'ylabel'),'units','normalized','position',[-0.19 0.53],'Fontsize',FntSz+6)
text(-0.175,0.555,'(T_2)','Fontsize',FntSz+2,'units','normalized')
text(-0.2,0.4,'(s^{-1/2})','Fontsize',FntSz,'units','normalized')
text(-0.25,1.05,'d','Fontsize',FntSz+8,'units','normalized','Fontweight','bold')
text(1.66,1.1,'                ','Fontsize',FntSz,...
    'HorizontalAlignment','center','Edgecolor',[0 0 0])
text(1.33,1.11,'\color{red}+','Fontsize',FntSz,...
    'HorizontalAlignment','right')
text(1.3,1.1,'\color{black} Outlier','Fontsize',FntSz,...
    'HorizontalAlignment','left')


% print('figure3','-dtiff','-r600')

