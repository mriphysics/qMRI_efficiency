%% Plot optimal results for the simple MRF experiments
%David Leitao (david.leitao@kcl.ac.uk); 23-04-20

clearvars; close all; clc;

%% Load optimal acquisition settings

N = 5;
T1 = 781;       T2 = 65;
dT1 = T1*1e-3;  dT2 = T2*1e-3;
TEmin = 2; TE0 = TEmin*ones(N,1); RF0 = zeros(N,1);

% non-DE results
load('simple_experiment_opt_acq_set_nonDE_N5.mat')
nonDE_T1.Acq    = T1_estimation.acqSet{1};
nonDE_T1.Mxy    = T1_estimation.allSig{1};
nonDE_T2.Acq    = T2_estimation.acqSet{1};
nonDE_T2.Mxy    = T2_estimation.allSig{1};
nonDE_joint.Acq = joint_estimation.acqSet{1};
nonDE_joint.Mxy = joint_estimation.allSig{1};
% DE results
load('simple_experiment_opt_acq_set_DE_N5.mat')
DE_joint.Acq = joint_estimation.acqSet{1};
DE_joint.Mxy = joint_estimation.allSig{1};

% calculate signals and its derivatives for the optimal acquisition
% settings
% T1 scenario: dmdT1
[~, aux_signal1] = ...
        cppAnalysis_acq_set_spoiledMRF(N, nonDE_T1.Acq(1:N), RF0, nonDE_T1.Acq(N+1:2*N), TE0, T1+dT1, T2);
[~, aux_signal2] = ...
        cppAnalysis_acq_set_spoiledMRF(N, nonDE_T1.Acq(1:N), RF0, nonDE_T1.Acq(N+1:2*N), TE0, T1-dT1, T2);
nonDE_T1.dMxydT1 = (aux_signal1-aux_signal2) / (2*dT1);

% T2 scenario: dmdT2
[~, aux_signal1] = ...
        cppAnalysis_acq_set_spoiledMRF(N, nonDE_T2.Acq(1:N), RF0, nonDE_T2.Acq(N+1:2*N), TE0, T1, T2+dT2);
[~, aux_signal2] = ...
        cppAnalysis_acq_set_spoiledMRF(N, nonDE_T2.Acq(1:N), RF0, nonDE_T2.Acq(N+1:2*N), TE0, T1, T2-dT2);
nonDE_T2.dMxydT2 = (aux_signal1-aux_signal2) / (2*dT2);
   
% non-DE joint scenario: dmdT1
[~, aux_signal1] = ...
        cppAnalysis_acq_set_spoiledMRF(N, nonDE_joint.Acq(1:N), RF0, nonDE_joint.Acq(N+1:2*N), TE0, T1+dT1, T2);
[~, aux_signal2] = ...
        cppAnalysis_acq_set_spoiledMRF(N, nonDE_joint.Acq(1:N), RF0, nonDE_joint.Acq(N+1:2*N), TE0, T1-dT1, T2);
nonDE_joint.dMxydT1 = (aux_signal1-aux_signal2) / (2*dT1);
% non-DE joint scenario: dmdT2
[~, aux_signal1] = ...
        cppAnalysis_acq_set_spoiledMRF(N, nonDE_joint.Acq(1:N), RF0, nonDE_joint.Acq(N+1:2*N), TE0, T1, T2+dT2);
[~, aux_signal2] = ...
        cppAnalysis_acq_set_spoiledMRF(N, nonDE_joint.Acq(1:N), RF0, nonDE_joint.Acq(N+1:2*N), TE0, T1, T2-dT2);
nonDE_joint.dMxydT2 = (aux_signal1-aux_signal2) / (2*dT2);

% DE joint scenario: dmdT1
[~, aux_signal1] = ...
        cppAnalysis_acq_set_DE_spoiledMRF(N, DE_joint.Acq(1:N), RF0, DE_joint.Acq(N+1:2*N), TE0, T1+dT1, T2);
[~, aux_signal2] = ...
        cppAnalysis_acq_set_DE_spoiledMRF(N, DE_joint.Acq(1:N), RF0, DE_joint.Acq(N+1:2*N), TE0, T1-dT1, T2);
DE_joint.dMxydT1 = (aux_signal1-aux_signal2) / (2*dT1);
% DE joint scenario: dmdT2
[~, aux_signal1] = ...
        cppAnalysis_acq_set_DE_spoiledMRF(N, DE_joint.Acq(1:N), RF0, DE_joint.Acq(N+1:2*N), TE0, T1, T2+dT2);
[~, aux_signal2] = ...
        cppAnalysis_acq_set_DE_spoiledMRF(N, DE_joint.Acq(1:N), RF0, DE_joint.Acq(N+1:2*N), TE0, T1, T2-dT2);
DE_joint.dMxydT2 = (aux_signal1-aux_signal2) / (2*dT2);    
    
clear T1_estimation T2_estimation joint_estimation aux_signal1 aux_signal2

%% Plot results for simple MRF experiments

figure;
set(gcf,'Units','normalized','Outerposition',[0 0 1 1],'Color','w')

dr     = 0.03; yarrow = 0.70;
LW = 2; FntSz = 16; freq = 1.5;
f_sinc =@(npts) 0.2 * sin(freq*linspace(-20,20,npts)) .* hanning(npts)';
cgray = 0.55 * [1 1 1]; c = lines(7);
TE_x = 0.05; TE = 0;
blueArrow = [0 0.6 0.2];

% subplot for sequence structure 
h1 = subplot(5,1,1);
hold on; 
% zero RF
plot([0.5 0.5], [0 yarrow], 'Linewidth', LW+1, 'Color', cgray)
plot([-dr 0 dr]+0.5, [yarrow-2*dr yarrow yarrow-2*dr], 'Linewidth', LW+1, 'Color', cgray)
text(0.5, yarrow+12*dr, {['RF pulse']; ['(\alpha_{n-1}, \phi_{n-1})']}, 'Fontsize', FntSz, 'HorizontalAlignment', 'center', 'Color', cgray)
% zero gradient
area([1.45 1.5 1.8 1.85],[0 0.3 0.3 0], 'FaceColor', cgray, 'EdgeColor', cgray)
text(1.65,0.5,'G_{spoiling}', 'Fontsize', FntSz, 'HorizontalAlignment', 'center', 'Color', cgray)
% first RF
plot([2 2], [0 yarrow], '-k', 'Linewidth', LW+1)
plot([2-dr 2 2+dr], [yarrow-2*dr yarrow yarrow-2*dr], '-k', 'Linewidth', LW+1)
text(2, yarrow+12*dr, {['RF pulse']; ['(\alpha_n, \phi_n)']}, 'Fontsize', FntSz, 'HorizontalAlignment', 'center')
% first gradient
area([2.95 3 3.3 3.35],[0 0.3 0.3 0], 'FaceColor', [0 0 0], 'EdgeColor', [0 0 0])
text(3.15,0.5,'G_{spoiling}', 'Fontsize', FntSz, 'HorizontalAlignment', 'center')

% second RF
plot([3.5 3.5], [0 yarrow], 'Linewidth', LW+1, 'Color', cgray)
plot([3.5-dr 3.5 3.5+dr], [yarrow-2*dr yarrow yarrow-2*dr], 'Linewidth', LW+1, 'Color', cgray)
text(3.5, yarrow+12*dr, {['RF pulse']; ['(\alpha_{n+1}, \phi_{n+1})']}, 'Fontsize', FntSz, 'HorizontalAlignment', 'center', 'Color', cgray)
% second gradient
area([4.45 4.5 4.8 4.85],[0 0.3 0.3 0], 'FaceColor', cgray, 'EdgeColor', cgray)
text(4.65,0.5,'G_{spoiling}', 'Fontsize', FntSz, 'HorizontalAlignment', 'center', 'Color', cgray)
%%%
plot([-0.5 6],[0 0], '-k', 'Linewidth', LW)
plot(6+[-dr 0 -dr],[-dr 0 dr], '-k', 'Linewidth', LW)
%%%
plot([-0.5 2],[0 0], '-', 'Color', cgray, 'Linewidth', LW)
plot([2    3.5],[0 0], '-k', 'Linewidth', LW)
plot([3.5    6],[0 0], '-', 'Color', cgray, 'Linewidth', LW)
plot(6+[-dr 0 -dr],[-dr 0 dr], '-', 'Color', cgray, 'Linewidth', LW)
%%%
text(6+dr,-2*dr, 't', 'FontWeight', 'bold', 'FontSize', FntSz)
% zero readout
plot(linspace(0.5+TE_x,0.9+TE_x,2000),f_sinc(2000), 'Color', cgray, 'Linewidth', LW-1)
text(0.7+TE_x,0.375,'readout', 'Fontsize', FntSz, 'HorizontalAlignment', 'center', 'Color', cgray)
% first readout
plot(linspace(2+TE_x,2.4+TE_x,2000),f_sinc(2000),'k', 'Linewidth', LW-1)
text(2.25+dr,0.375,'readout', 'Fontsize', FntSz, 'HorizontalAlignment', 'center')
% second readout
plot(linspace(3.5+TE_x,3.9+TE_x,2000),f_sinc(2000), 'Color', cgray, 'Linewidth', LW-1)
text(3.7+dr,0.375,'readout', 'Fontsize', FntSz, 'HorizontalAlignment', 'center', 'Color', cgray)
%%% 3 dots
plot([-0.05 0 0.05],yarrow/2*ones(3,1),'.', 'Color', cgray, 'Linewidth', LW+1)
plot([5.25 5.30 5.35],yarrow/2*ones(3,1),'.', 'Color', cgray, 'Linewidth', LW+1)
axis([-0.5 6 -0.2 1])
set(gca,'visible','off')
text(-0.05,1.15,'a','Units','normalized','Fontsize',20,'FontWeight','bold')



% subplot for T1
h2 = subplot(5,1,2);
hold on; 
plot([-0.5 5.75],[0 0], '-k', 'Linewidth', LW)
plot(5.75+[-dr 0 -dr],[-dr 0 dr], '-k', 'Linewidth', LW)
text(5.75+dr,-2*dr, 't (ms)', 'FontWeight', 'bold', 'FontSize', FntSz)
xpos = [0.5 2.3 3.7 4.35 5];
plot(linspace(-0.495,5.5,round(5.5*20/yarrow)),yarrow*ones(round(5.5*20/yarrow),1),'k-', 'Linewidth', LW-1, 'Color',1.3*cgray)
plot(linspace(-0.495,5.5,round(5.5*20/yarrow)),yarrow/2*ones(round(5.5*20/yarrow),1),'-', 'Linewidth', 0.5, 'Color', 1.3*cgray)
text(-0.5-0.07,yarrow,'1', 'Fontsize', 16, 'HorizontalAlignment','right', 'Color',0.8*cgray)
plot([-0.5 -0.5+dr], [yarrow yarrow]./2, '-k', 'Linewidth', LW, 'Color',1.3*cgray)
text(-0.5-0.07,yarrow/2,'0.5', 'Fontsize', 16, 'HorizontalAlignment','right', 'Color',0.8*cgray)
text(-0.5-0.07,0,'0', 'Fontsize', 16, 'HorizontalAlignment','right', 'Color',0.8*cgray)
plot([-0.5 -0.5], [0 yarrow], '-k', 'Linewidth', LW, 'Color',1.3*cgray)
for ii=1:N
    sf = nonDE_T1.Acq(ii)/pi;
    plot(xpos(ii)*ones(20,1), linspace(0,yarrow,20),'.', 'Linewidth', 0.5, 'Color', cgray)
    if sf>0
        plot([xpos(ii) xpos(ii)], sf*[0 yarrow], '-', 'Linewidth', LW+1, 'Color', blueArrow)
        plot([xpos(ii)-dr xpos(ii) xpos(ii)+dr], [sf*yarrow-2*dr sf*yarrow sf*yarrow-2*dr], '-', 'Linewidth', LW+1, 'Color', blueArrow)
    end
    text(xpos(ii), yarrow+5*dr, ['\alpha_',num2str(ii),'=',num2str(rad2deg(nonDE_T1.Acq(ii)),'%.0f'),'\circ'], ...
        'Fontsize', 16, 'HorizontalAlignment','center')
    if ii==1
        text(xpos(ii)+TE_x, -8*dr, ['t_',num2str(ii),'=',num2str(TE)], 'Fontsize', 16, 'HorizontalAlignment','center')
    else
        text(xpos(ii)+TE_x, -8*dr, ['t_',num2str(ii),'=',num2str(sum(nonDE_T1.Acq(N+1:N+ii-1)+TE),'%.0f')], 'Fontsize', 16, 'HorizontalAlignment','center')
    end
end
plot(xpos+TE_x, yarrow*abs(nonDE_T1.Mxy), 'o', 'Linewidth', LW-1, 'MarkerFaceColor',c(1,:), 'MarkerEdgeColor', c(1,:), 'Color', c(1,:))
plot(xpos+TE_x, yarrow*T1*abs(nonDE_T1.dMxydT1), 's', 'Linewidth', LW-1, 'MarkerFaceColor',c(2,:), 'MarkerEdgeColor', c(2,:), 'Color', c(2,:))
axis([-0.5 5.75 -dr 1])
set(gca,'visible','off')
set(gca, 'Fontsize', FntSz)
hlines = findobj(h2,'Type', 'line');
hlegend = legend([hlines(2:-1:1)],'s','\partials/\partialT1',...
    'Orientation','Vertical','Location','best');
hlegend.Position(1) = 0.85;
hlegend.Position(2) = hlegend.Position(2)-0.025; 
text(0,yarrow/2,'M(t_0)=M_0','Fontsize',16,'HorizontalAlignment','center','Color',0.8*cgray,'BackgroundColor','w')
text(-0.05,1.1,'b','Units','normalized','Fontsize',20,'FontWeight','bold')



% subplot for T2
h3 = subplot(5,1,3);
hold on; 
plot([-0.5 5.75],[0 0], '-k', 'Linewidth', LW)
plot(5.75+[-dr 0 -dr],[-dr 0 dr], '-k', 'Linewidth', LW)
text(5.75+dr,-2*dr, 't (ms)', 'FontWeight', 'bold', 'FontSize', FntSz)
xpos = [0.5 1.4 3.2 4.1 5];
plot(linspace(-0.495,5.5,round(5.5*20/yarrow)),yarrow*ones(round(5.5*20/yarrow),1),'k-', 'Linewidth', LW-1, 'Color',1.3*cgray)
plot(linspace(-0.495,5.5,round(5.5*20/yarrow)),yarrow/2*ones(round(5.5*20/yarrow),1),'-', 'Linewidth', 0.5, 'Color', 1.3*cgray)
text(-0.5-0.07,yarrow,'1', 'Fontsize', 16, 'HorizontalAlignment','right', 'Color',0.8*cgray)
plot([-0.5 -0.5+dr], [yarrow yarrow]./2, '-k', 'Linewidth', LW, 'Color',1.3*cgray)
text(-0.5-0.07,yarrow/2,'0.5', 'Fontsize', 16, 'HorizontalAlignment','right', 'Color',0.8*cgray)
text(-0.5-0.07,0,'0', 'Fontsize', 16, 'HorizontalAlignment','right', 'Color',0.8*cgray)
plot([-0.5 -0.5], [0 yarrow], '-k', 'Linewidth', LW, 'Color',1.3*cgray)
for ii=1:N
    sf = nonDE_T2.Acq(ii)/pi;
    plot(xpos(ii)*ones(20,1), linspace(0,yarrow,20),'.', 'Linewidth', 0.5, 'Color', cgray)
    if sf>0
        plot([xpos(ii) xpos(ii)], sf*[0 yarrow], '-', 'Linewidth', LW+1, 'Color', blueArrow)
        plot([xpos(ii)-dr xpos(ii) xpos(ii)+dr], [sf*yarrow-2*dr sf*yarrow sf*yarrow-2*dr], '-', 'Linewidth', LW+1, 'Color', blueArrow)
    end
    text(xpos(ii), yarrow+5*dr, ['\alpha_',num2str(ii),'=',num2str(rad2deg(nonDE_T2.Acq(ii)),'%.0f'),'\circ'], ...
        'Fontsize', 16, 'HorizontalAlignment','center')
    if ii==1
        text(xpos(ii)+TE_x, -8*dr, ['t_',num2str(ii),'=',num2str(TE)], 'Fontsize', 16, 'HorizontalAlignment','center')
    else
        text(xpos(ii)+TE_x, -8*dr, ['t_',num2str(ii),'=',num2str(sum(nonDE_T2.Acq(N+1:N+ii-1)+TE),'%.0f')], 'Fontsize', 16, 'HorizontalAlignment','center')
    end
end
plot(xpos+TE_x, yarrow*abs(nonDE_T2.Mxy), 'o', 'Linewidth', LW-1, 'MarkerFaceColor',c(1,:), 'MarkerEdgeColor', c(1,:), 'Color', c(1,:))
plot(xpos+TE_x, yarrow*T2*abs(nonDE_T2.dMxydT2), 'd', 'Linewidth', LW-1, 'MarkerFaceColor',c(3,:), 'MarkerEdgeColor', c(3,:), 'Color', c(3,:))
axis([-0.5 5.75 -dr 1])
set(gca,'visible','off')
set(gca, 'Fontsize', FntSz)
hlines = findobj(h3,'Type', 'line');
hlegend = legend([hlines(2:-1:1)],'s','\partials/\partialT2',...
    'Orientation','Vertical','Location','best');
hlegend.Position(1) = 0.85;
hlegend.Position(2) = hlegend.Position(2)-0.025; 
text(0,yarrow/2,'M(t_0)=M_0','Fontsize',16,'HorizontalAlignment','center','Color',0.8*cgray,'BackgroundColor','w')
text(-0.05,1.1,'c','Units','normalized','Fontsize',20,'FontWeight','bold')



% subplot for non-DE joint
h4 = subplot(5,1,4);
hold on; 
plot([-0.5 5.75],[0 0], '-k', 'Linewidth', LW)
plot(5.75+[-dr 0 -dr],[-dr 0 dr], '-k', 'Linewidth', LW)
text(5.75+dr,-2*dr, 't (ms)', 'FontWeight', 'bold', 'FontSize', FntSz)
xpos = [0.5 1.23 2.75 3.875 5]; 
plot(linspace(-0.495,5.5,round(5.5*20/yarrow)),yarrow*ones(round(5.5*20/yarrow),1),'k-', 'Linewidth', LW-1, 'Color',1.3*cgray)
plot(linspace(-0.495,5.5,round(5.5*20/yarrow)),yarrow/2*ones(round(5.5*20/yarrow),1),'-', 'Linewidth', 0.5, 'Color', 1.3*cgray)
text(-0.5-0.07,yarrow,'1', 'Fontsize', 16, 'HorizontalAlignment','right', 'Color',0.8*cgray)
plot([-0.5 -0.5+dr], [yarrow yarrow]./2, '-k', 'Linewidth', LW, 'Color',1.3*cgray)
text(-0.5-0.07,yarrow/2,'0.5', 'Fontsize', 16, 'HorizontalAlignment','right', 'Color',0.8*cgray)
text(-0.5-0.07,0,'0', 'Fontsize', 16, 'HorizontalAlignment','right', 'Color',0.8*cgray)
plot([-0.5 -0.5], [0 yarrow], '-k', 'Linewidth', LW, 'Color',1.3*cgray)
for ii=1:N
    sf = nonDE_joint.Acq(ii)/pi;
    plot(xpos(ii)*ones(20,1), linspace(0,yarrow,20),'.', 'Linewidth', 0.5, 'Color', cgray)
    if sf>0
        plot([xpos(ii) xpos(ii)], sf*[0 yarrow], '-', 'Linewidth', LW+1, 'Color', blueArrow)
        plot([xpos(ii)-dr xpos(ii) xpos(ii)+dr], [sf*yarrow-2*dr sf*yarrow sf*yarrow-2*dr], '-', 'Linewidth', LW+1, 'Color', blueArrow)
    end
    text(xpos(ii), yarrow+5*dr, ['\alpha_',num2str(ii),'=',num2str(rad2deg(nonDE_joint.Acq(ii)),'%.0f'),'\circ'], ...
        'Fontsize', 16, 'HorizontalAlignment','center')
    if ii==1
        text(xpos(ii)+TE_x, -8*dr, ['t_',num2str(ii),'=',num2str(TE)], 'Fontsize', 16, 'HorizontalAlignment','center')
    else
        text(xpos(ii)+TE_x, -8*dr, ['t_',num2str(ii),'=',num2str(sum(nonDE_joint.Acq(N+1:N+ii-1)+TE),'%.0f')], 'Fontsize', 16, 'HorizontalAlignment','center')
    end
end
plot(xpos+TE_x, yarrow*abs(nonDE_joint.Mxy), 'o', 'Linewidth', LW-1, 'MarkerFaceColor',c(1,:), 'MarkerEdgeColor', c(1,:), 'Color', c(1,:))
plot(xpos+TE_x, yarrow*T1*abs(nonDE_joint.dMxydT1), 's', 'Linewidth', LW-1, 'MarkerFaceColor',c(2,:), 'MarkerEdgeColor', c(2,:), 'Color', c(2,:))
plot(xpos+TE_x, yarrow*T2*abs(nonDE_joint.dMxydT2), 'd', 'Linewidth', LW-1, 'MarkerFaceColor',c(3,:), 'MarkerEdgeColor', c(3,:), 'Color', c(3,:))
axis([-0.5 5.75 -dr 1])
set(gca,'visible','off')
set(gca, 'Fontsize', FntSz)
hlines = findobj(h4,'Type', 'line');
hlegend = legend([hlines(3:-1:1)],'s','\partials/\partialT1','\partials/\partialT2',...
    'Orientation','Vertical','Location','best');
hlegend.Position(1) = 0.85;
hlegend.Position(2) = hlegend.Position(2)-0.025; 
text(0,yarrow/2,'M(t_0)=M_0','Fontsize',16,'HorizontalAlignment','center','Color',0.8*cgray,'BackgroundColor','w')
text(-0.05,1.1,'d','Units','normalized','Fontsize',20,'FontWeight','bold')



% subplot for DE joint
h5 = subplot(5,1,5);
hold on; 
plot([-0.5 5.75],[0 0], '-k', 'Linewidth', LW)
plot(5.75+[-dr 0 -dr],[-dr 0 dr], '-k', 'Linewidth', LW)
text(5.75+dr,-2*dr, 't (ms)', 'FontWeight', 'bold', 'FontSize', FntSz)
xpos = [0.5 1 2 2.5 4.2]; 
plot(linspace(-0.495,5.5,round(5.5*20/yarrow)),yarrow*ones(round(5.5*20/yarrow),1),'k-', 'Linewidth', LW-1, 'Color',1.3*cgray)
plot(linspace(-0.495,5.5,round(5.5*20/yarrow)),yarrow/2*ones(round(5.5*20/yarrow),1),'-', 'Linewidth', 0.5, 'Color', 1.3*cgray)
text(-0.5-0.07,yarrow,'1', 'Fontsize', 16, 'HorizontalAlignment','right', 'Color',0.8*cgray)
plot([-0.5 -0.5+dr], [yarrow yarrow]./2, '-k', 'Linewidth', LW, 'Color',1.3*cgray)
text(-0.5-0.07,yarrow/2,'0.5', 'Fontsize', 16, 'HorizontalAlignment','right', 'Color',0.8*cgray)
text(-0.5-0.07,0,'0', 'Fontsize', 16, 'HorizontalAlignment','right', 'Color',0.8*cgray)
plot([-0.5 -0.5], [0 yarrow], '-k', 'Linewidth', LW, 'Color',1.3*cgray)
text(0.1,0,'//', 'Fontsize', 16, 'HorizontalAlignment','center')
text(4.6,0,'//', 'Fontsize', 16, 'HorizontalAlignment','center')
for ii=1:N
    sf = DE_joint.Acq(ii)/pi;
    plot(xpos(ii)*ones(20,1), linspace(0,yarrow,20),'.', 'Linewidth', 0.5, 'Color', cgray)
    if sf>0
        plot([xpos(ii) xpos(ii)], sf*[0 yarrow], '-', 'Linewidth', LW+1, 'Color', blueArrow)
        plot([xpos(ii)-dr xpos(ii) xpos(ii)+dr], [sf*yarrow-2*dr sf*yarrow sf*yarrow-2*dr], '-', 'Linewidth', LW+1, 'Color', blueArrow)
    end
    text(xpos(ii), yarrow+5*dr, ['\alpha_',num2str(ii),'=',num2str(rad2deg(DE_joint.Acq(ii)),'%.0f'),'\circ'], ...
        'Fontsize', 16, 'HorizontalAlignment','center')
    if ii==1
        text(xpos(ii)+TE_x, -8*dr, ['t_',num2str(ii),'=',num2str(TE)], 'Fontsize', 16, 'HorizontalAlignment','center')
    else
        text(xpos(ii)+TE_x, -8*dr, ['t_',num2str(ii),'=',num2str(sum(DE_joint.Acq(N+1:N+ii-1)+TE),'%.0f')], 'Fontsize', 16, 'HorizontalAlignment','center')
    end
end
plot(5*ones(20,1), linspace(0,yarrow,20),'.', 'Linewidth', 0.5, 'Color', cgray)
text(5, yarrow+5*dr, ['\alpha_',num2str(1),'=',num2str(rad2deg(DE_joint.Acq(1)),'%.0f'),'\circ'], ...
        'Fontsize', 16, 'HorizontalAlignment','center', 'Color', cgray)    
sf = DE_joint.Acq(1)/pi;
plot([5 5], sf*[0 yarrow], '-', 'Linewidth', LW+1, 'Color', mat2gray(blueArrow+0.4))
plot([5-dr 5 5+dr], [sf*yarrow-2*dr sf*yarrow sf*yarrow-2*dr], '-', 'Linewidth', LW+1, 'Color', mat2gray(blueArrow+0.4))
text(5+TE_x, -8*dr, ['t_',num2str(1),'=',num2str(sum(DE_joint.Acq(N+1:N+5)+TE),'%.0f')], 'Fontsize', 16, 'HorizontalAlignment','center', 'Color', cgray)
plot(5+TE_x, yarrow*abs(DE_joint.Mxy(1)), 'o', 'Linewidth', LW-1, 'MarkerFaceColor',mat2gray(c(1,:)+0.3), 'MarkerEdgeColor', mat2gray(c(1,:)+0.4), 'Color', c(1,:))
plot(5+TE_x, yarrow*T1*abs(DE_joint.dMxydT1(1)), 's', 'Linewidth', LW-1, 'MarkerFaceColor',mat2gray(c(2,:)+0.3), 'MarkerEdgeColor', mat2gray(c(2,:)+0.4), 'Color', c(2,:))
plot(5+TE_x, yarrow*T2*abs(DE_joint.dMxydT2(1)), 'd', 'Linewidth', LW-1, 'MarkerFaceColor',mat2gray(c(3,:)+0.3), 'MarkerEdgeColor', mat2gray(c(3,:)+0.4), 'Color', c(3,:))
plot(xpos+TE_x, yarrow*abs(DE_joint.Mxy), 'o', 'Linewidth', LW-1, 'MarkerFaceColor',c(1,:), 'MarkerEdgeColor', c(1,:), 'Color', c(1,:))
plot(xpos+TE_x, yarrow*T1*abs(DE_joint.dMxydT1), 's', 'Linewidth', LW-1, 'MarkerFaceColor',c(2,:), 'MarkerEdgeColor', c(2,:), 'Color', c(2,:))
plot(xpos+TE_x, yarrow*T2*abs(DE_joint.dMxydT2), 'd', 'Linewidth', LW-1, 'MarkerFaceColor',c(3,:), 'MarkerEdgeColor', c(3,:), 'Color', c(3,:))
axis([-0.5 5.75 -dr 1])
set(gca,'visible','off')
set(gca, 'Fontsize', FntSz)
hlines = findobj(h5,'Type', 'line');
hlegend = legend([hlines(3:-1:1)],'s','\partials/\partialT1','\partials/\partialT2',...
    'Orientation','Vertical','Location','best');
hlegend.Position(1) = 0.85;
hlegend.Position(2) = hlegend.Position(2)-0.025; 
text(-0.05,1.1,'e','Units','normalized','Fontsize',20,'FontWeight','bold')


% print('figure2','-dtiff','-r600')

