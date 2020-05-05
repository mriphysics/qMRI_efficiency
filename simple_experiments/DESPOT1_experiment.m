%% Optimisation of DESPOT1
%David Leitao (david.leitao@kcl.ac.uk); 23-04-20

%Requires the handles to the SPGR steady-state signal function and its 
%derivative -> needs path to library folder

clearvars; close all; clc;

%% Sequence constraints and tissue parameters values

TRmin =  5; %Minimum TR [ms]
TEspgr = 2; %SPGR echo time [ms]

T1 = 781;   %Spin-lattice relaxation constant [ms]
T2 = 65;    %Spin-spin relaxation constant [ms] 
M0 = 1;     %Equilibrium magnetisation [a.u.]
P0 = 0;     %Signal phase [rad]
B0 = 0;     %Off-resonance [Hz]

%% Loads signal function and derivatives 
% vector u is the acquisition settings and vector p is the tissue parameters
%   u = [FA, TR]
%   p = [T1, T2, M0]

% load handles to calculate steady-state signals
load('hnd_SignalFunctions_SteadyState.mat')
signal_SPGR  =@(u,p) func_signal_SPGR(B0,u(1),p(3),P0,p(1),p(2),TEspgr,u(2));

% load handles to calculate signal derivatives
load('hnd_SignalDerivatives_SteadyState.mat')
dmdT1_SPGR  =@(u,p) func_dmdT1_SPGR(B0,u(1),p(3),P0,p(1),p(2),TEspgr,u(2)); 
dmdM0_SPGR  =@(u,p) func_dmdM0_SPGR(B0,u(1),P0,p(1),p(2),TEspgr,u(2));

%% Optimisation options

options = optimoptions('fmincon',...
    'SpecifyConstraintGradient',false,...
    'Algorithm','sqp',...
    'MaxFunctionEvaluations',Inf,...
    'Display','none',...
    'OptimalityTolerance',1e-4,...
    'StepTolerance',1e-4,...
    'MaxIterations',1e4,...
    'Useparallel',false);

%% Run optimisation

% define number of SPGR signals
nSPGR  = 2; 

func =@(u,p) cost_function_DESPOT1(u, p, nSPGR, dmdT1_SPGR,  dmdM0_SPGR);
CostFunc =@(x) func([x(1) x(2) x(3) x(3)], [T1 T2 M0]); %Same TR for both SPGR

% constraints for current combination
lb = [0*ones(nSPGR,1);    TRmin];
ub = [pi/2*ones(nSPGR,1); Inf];
% bounds used to create random initialisations (Inf cannot be used)
lb0 = [0*ones(nSPGR,1);    TRmin];
ub0 = [pi/2*ones(nSPGR,1); 100];

nMS = 100;
rng('default')

all_fval = Inf*ones(nMS, 1);
all_uopt = zeros(nMS, length(ub));

tic
for ii=1:nMS
    % random initialisation within the previously defined bounds   
    u0 = 0.8*rand(nSPGR+1,1).*(ub0-lb0) + 1.1*lb0;
    
    % try random strating point; if it fails, try next one
    try
        [all_uopt(ii,:), all_fval(ii)] = fmincon(CostFunc, u0, [], [], [], [], lb, ub, [], options);
    catch
        continue;
    end
end
toc
% extract best solution
all_fval(all_fval<0) = Inf;
idx_best = find(min(all_fval)==all_fval,1,'first');
uopt = all_uopt(idx_best,:);

eff_T1_DESPOT1 = 1 / sqrt(CostFunc(uopt));
fprintf('\nDESPOT1 T1 efficiency = %.4f \n\n',eff_T1_DESPOT1)
fprintf('%c1 = %.1f deg \t %c1 = %.1f deg \n',char(945),rad2deg(uopt(1)),char(945),rad2deg(uopt(2)))
fprintf('TR = %.1f ms \n\n',uopt(3))

%% Plot heat maps of T1-to-noise ratio, T1 efficiency and optimal acq. set.

c = [0 0.7 0];

figure; 
set(gcf,'Units','normalized','outerposition',[0 0.25 1 0.45],'Color','w')

max_y = 0.06;
FntSz = 18;

minFA  = 0;
maxFA  = 20; 
FAstep = 0.1;
minTR  = 2;
maxTR  = 20;
TRstep = 0.2;

optFA1 = rad2deg(uopt(1));
optFA2 = rad2deg(uopt(2));

cutA1 = floor((1/FAstep)*(optFA1))*FAstep;
cutA2 = floor((1/FAstep)*(optFA2))*FAstep;
cutA1 = min([cutA1 cutA2]);
cutA2 = cutA1;

cutTR = floor((1/TRstep)*(uopt(3)))*TRstep;


SNRmax = 1000;

hsp2 = subplot(1,3,1);
text(-0.25,1.05,'a','Units','normalized','Fontsize',FntSz+6,'FontWeight','bold')


%%%%%%%%%%%%%%%%%%%%%%%%%%%% A1-TR plane
[auxTR, auxA1, auxA2] = meshgrid(cutTR:TRstep:maxTR, cutA1:FAstep:maxFA, cutA2);
auxC = zeros(size(auxTR));
for ii=1:numel(auxTR)
    auxC(ii) = SNRmax*sqrt(auxTR(ii)*1e-3)/sqrt(CostFunc([auxA1(ii)*pi/180, auxA2(ii)*pi/180, auxTR(ii), auxTR(ii)])); 
end

hold on; 
s = surface(squeeze(auxA1), squeeze(auxA2), squeeze(auxTR), squeeze(auxC));
s.EdgeColor = 'none';
s.FaceColor = 'interp';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%% A2-TR plane
[auxTR, auxA1, auxA2] = meshgrid(cutTR:TRstep:maxTR, maxFA, minFA:FAstep:cutA2);
auxC = zeros(size(auxTR));
for ii=1:numel(auxTR)
    auxC(ii) = SNRmax*sqrt(auxTR(ii)*1e-3)/sqrt(CostFunc([auxA1(ii)*pi/180, auxA2(ii)*pi/180, auxTR(ii), auxTR(ii)])); 
end
s = surface(squeeze(auxA1), squeeze(auxA2), squeeze(auxTR), squeeze(auxC));
s.EdgeColor = 'none';
s.FaceColor = 'interp';

[auxTR, auxA1, auxA2] = meshgrid(minTR:TRstep:cutTR, maxFA, minFA:FAstep:maxFA);
auxC = zeros(size(auxTR));
for ii=1:numel(auxTR)
    auxC(ii) = SNRmax*sqrt(auxTR(ii)*1e-3)/sqrt(CostFunc([auxA1(ii)*pi/180, auxA2(ii)*pi/180, auxTR(ii), auxTR(ii)])); 
end
s = surface(squeeze(auxA1), squeeze(auxA2), squeeze(auxTR), squeeze(auxC));
s.EdgeColor = 'none';
s.FaceColor = 'interp';
% s.FaceAlpha = Falpha;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%% A1-A2 plane
[auxTR, auxA1, auxA2] = meshgrid(cutTR, cutA1:FAstep:maxFA, cutA2:FAstep:maxFA);
auxC = zeros(size(auxTR));
for ii=1:numel(auxTR)
    if auxA1(ii)+FAstep>=auxA2(ii)
        auxC(ii) = SNRmax*sqrt(auxTR(ii)*1e-3)/sqrt(CostFunc([auxA1(ii)*pi/180, auxA2(ii)*pi/180, auxTR(ii), auxTR(ii)])); 
    else
        auxC(ii) = NaN;
    end
end

s = surface(squeeze(auxA1), squeeze(auxA2), squeeze(auxTR), squeeze(auxC));
s.EdgeColor = 'none';
s.FaceColor = 'interp';

[auxTR, auxA1, auxA2] = meshgrid(maxTR, minFA:FAstep:cutA1, minFA:FAstep:maxFA);
auxC = zeros(size(auxTR));
for ii=1:numel(auxTR)
    if auxA1(ii)+FAstep>=auxA2(ii)
        auxC(ii) = SNRmax*sqrt(auxTR(ii)*1e-3)/sqrt(CostFunc([auxA1(ii)*pi/180, auxA2(ii)*pi/180, auxTR(ii), auxTR(ii)])); 
    else
        auxC(ii) = NaN;
    end
end

s = surface(squeeze(auxA1), squeeze(auxA2), squeeze(auxTR), squeeze(auxC));
s.EdgeColor = 'none';
s.FaceColor = 'interp';

[auxTR, auxA1, auxA2] = meshgrid(maxTR, cutA1:FAstep:maxFA, minFA:FAstep:cutA2);
auxC = zeros(size(auxTR));
for ii=1:numel(auxTR)
    auxC(ii) = SNRmax*sqrt(auxTR(ii)*1e-3)/sqrt(CostFunc([auxA1(ii)*pi/180, auxA2(ii)*pi/180, auxTR(ii), auxTR(ii)])); 
end

s = surface(squeeze(auxA1), squeeze(auxA2), squeeze(auxTR), squeeze(auxC));
s.EdgeColor = 'none';
s.FaceColor = 'interp';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%


plot3([minFA maxFA],[minFA minFA],[maxTR maxTR],'-k','Linewidth',1)
plot3([cutA1 maxFA],[cutA2 cutA2],[maxTR maxTR],'-k','Linewidth',1)
plot3([cutA1 maxFA],[cutA2 cutA2],[cutTR cutTR],'-k','Linewidth',1)
plot3([maxFA maxFA],[minFA cutA2],[maxTR maxTR],'-k','Linewidth',1)
plot3([maxFA maxFA],[minFA maxFA],[minTR minTR],'-k','Linewidth',1)
plot3([maxFA maxFA],[cutA2 maxFA],[cutTR cutTR],'-k','Linewidth',1)
plot3([cutA1 cutA1],[cutA2 cutA2],[cutTR maxTR],'-k','Linewidth',1)
plot3([maxFA maxFA],[cutA2 cutA2],[cutTR maxTR],'-k','Linewidth',1)
plot3([cutA1 maxFA],[cutA2 maxFA],[cutTR cutTR],'-k','Linewidth',1)
plot3([minFA cutA1],[minFA cutA2],[maxTR maxTR],'-k','Linewidth',1)
plot3([maxFA maxFA],[maxFA maxFA],[minTR cutTR],'-k','Linewidth',1)

plot3([minFA minFA],[minFA maxFA],[minTR minTR],'-k','Linewidth',0.5)
plot3([minFA 1],[15 15],[minTR minTR],'-k','Linewidth',0.5)
plot3([minFA 1],[10 10],[minTR minTR],'-k','Linewidth',0.5)
plot3([minFA 1],[5  5], [minTR minTR],'-k','Linewidth',0.5)

view([120 30])

colormap('inferno')
hcb2 = colorbar; caxis([5 20]);
title(hcb2,'T_1NR','Fontsize',FntSz+2)
hcb2.TickLabels{1} = '\leq5';

axis([minFA maxFA minFA maxFA minTR maxTR])
xticks([5 10 15 20]); xlabel('\alpha_2 (deg)')
yticks([5 10 15 20]); ylabel('\alpha_1 (deg)')
zticks([5 10 15 20]); zlabel('TR (ms)')

set(gca,'Fontsize',FntSz)


hcb2.Position(1) = hcb2.Position(1)-0.04;
hcb2.Position(2) = hcb2.Position(2)+0.15;
hcb2.Position(4) = hcb2.Position(4)-0.3;

hsp2.Position(1) = hsp2.Position(1)-0.06;
hsp2.Position(2) = hsp2.Position(2)+0.04;
hsp2.Position(3) = hsp2.Position(3)+0.02;
hsp2.Position(4) = hsp2.Position(4)-0.05;


set(gca,'SortMethod','ChildOrder');




minFA  = 0;
maxFA  = 20; 
FAstep = 0.1; 
minTR  = 2;
maxTR  = 20;
TRstep = 0.2; 

optFA1 = rad2deg(uopt(1));
optFA2 = rad2deg(uopt(2));

cutA1 = floor((1/FAstep)*(optFA1))*FAstep;
cutA2 = floor((1/FAstep)*(optFA2))*FAstep;
true_cutA1 = cutA1;
true_cutA2 = cutA2;
cutA1 = min([cutA1 cutA2]);
cutA2 = cutA1;

cutTR = floor((1/TRstep)*(uopt(3)))*TRstep;
Falpha = 0.5;

hsp1 = subplot(1,3,2);
text(-0.25,1.05,'b','Units','normalized','Fontsize',FntSz+6,'FontWeight','bold')

%%%%%%%%%%%%%%%%%%%%%%%%%%%% A1-TR plane
[auxTR, auxA1, auxA2] = meshgrid(cutTR:TRstep:maxTR, cutA1:FAstep:maxFA, cutA2);
auxC = zeros(size(auxTR));
for ii=1:numel(auxTR)
    auxC(ii) = 1/sqrt(CostFunc([auxA1(ii)*pi/180, auxA2(ii)*pi/180, auxTR(ii), auxTR(ii)])); 
end

hold on; 
s = surface(squeeze(auxA1), squeeze(auxA2), squeeze(auxTR), squeeze(auxC));
s.EdgeColor = 'none';
s.FaceColor = 'interp';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%% A2-TR plane
[auxTR, auxA1, auxA2] = meshgrid(cutTR:TRstep:maxTR, maxFA, minFA:FAstep:cutA2);
auxC = zeros(size(auxTR));
for ii=1:numel(auxTR)
    auxC(ii) = 1/sqrt(CostFunc([auxA1(ii)*pi/180, auxA2(ii)*pi/180, auxTR(ii), auxTR(ii)])); 
end
s = surface(squeeze(auxA1), squeeze(auxA2), squeeze(auxTR), squeeze(auxC));
s.EdgeColor = 'none';
s.FaceColor = 'interp';

[auxTR, auxA1, auxA2] = meshgrid(minTR:TRstep:cutTR, maxFA, minFA:FAstep:maxFA);
auxC = zeros(size(auxTR));
for ii=1:numel(auxTR)
    auxC(ii) = 1/sqrt(CostFunc([auxA1(ii)*pi/180, auxA2(ii)*pi/180, auxTR(ii), auxTR(ii)])); 
end
s = surface(squeeze(auxA1), squeeze(auxA2), squeeze(auxTR), squeeze(auxC));
s.EdgeColor = 'none';
s.FaceColor = 'interp';
% s.FaceAlpha = Falpha;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%% A1-A2 plane
[auxTR, auxA1, auxA2] = meshgrid(cutTR, cutA1:FAstep:maxFA, cutA2:FAstep:maxFA);
auxC = zeros(size(auxTR));
for ii=1:numel(auxTR)
    if auxA1(ii)+FAstep>=auxA2(ii)
        auxC(ii) = 1/sqrt(CostFunc([auxA1(ii)*pi/180, auxA2(ii)*pi/180, auxTR(ii), auxTR(ii)])); 
    else
        auxC(ii) = NaN;
    end
end

s = surface(squeeze(auxA1), squeeze(auxA2), squeeze(auxTR), squeeze(auxC));
s.EdgeColor = 'none';
s.FaceColor = 'interp';

[auxTR, auxA1, auxA2] = meshgrid(maxTR, minFA:FAstep:cutA1, minFA:FAstep:maxFA);
auxC = zeros(size(auxTR));
for ii=1:numel(auxTR)
    if auxA1(ii)+FAstep>=auxA2(ii)
        auxC(ii) = 1/sqrt(CostFunc([auxA1(ii)*pi/180, auxA2(ii)*pi/180, auxTR(ii), auxTR(ii)])); 
    else
        auxC(ii) = NaN;
    end
end

s = surface(squeeze(auxA1), squeeze(auxA2), squeeze(auxTR), squeeze(auxC));
s.EdgeColor = 'none';
s.FaceColor = 'interp';

[auxTR, auxA1, auxA2] = meshgrid(maxTR, cutA1:FAstep:maxFA, minFA:FAstep:cutA2);
auxC = zeros(size(auxTR));
for ii=1:numel(auxTR)
    auxC(ii) = 1/sqrt(CostFunc([auxA1(ii)*pi/180, auxA2(ii)*pi/180, auxTR(ii), auxTR(ii)])); 
end

s = surface(squeeze(auxA1), squeeze(auxA2), squeeze(auxTR), squeeze(auxC));
s.EdgeColor = 'none';
s.FaceColor = 'interp';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%


plot3([minFA maxFA],[minFA minFA],[maxTR maxTR],'-k','Linewidth',1)
plot3([cutA1 maxFA],[cutA2 cutA2],[maxTR maxTR],'-k','Linewidth',1)
plot3([cutA1 maxFA],[cutA2 cutA2],[cutTR cutTR],'-k','Linewidth',1)
plot3([maxFA maxFA],[minFA cutA2],[maxTR maxTR],'-k','Linewidth',1)
plot3([maxFA maxFA],[minFA maxFA],[minTR minTR],'-k','Linewidth',1)
plot3([maxFA maxFA],[cutA2 maxFA],[cutTR cutTR],'-k','Linewidth',1)
plot3([cutA1 cutA1],[cutA2 cutA2],[cutTR maxTR],'-k','Linewidth',1)
plot3([maxFA maxFA],[cutA2 cutA2],[cutTR maxTR],'-k','Linewidth',1)
plot3([cutA1 maxFA],[cutA2 maxFA],[cutTR cutTR],'-k','Linewidth',1)
plot3([minFA cutA1],[minFA cutA2],[maxTR maxTR],'-k','Linewidth',1)
plot3([maxFA maxFA],[maxFA maxFA],[minTR cutTR],'-k','Linewidth',1)

plot3([minFA minFA],[minFA maxFA],[minTR minTR],'-k','Linewidth',0.5)
plot3([minFA 1],[15 15],[minTR minTR],'-k','Linewidth',0.5)
plot3([minFA 1],[10 10],[minTR minTR],'-k','Linewidth',0.5)
plot3([minFA 1],[5  5], [minTR minTR],'-k','Linewidth',0.5)

view([120 30])

colormap('inferno')
hcb1 = colorbar; caxis([0.05 0.2]);
title(hcb1,'\eta(T_1)','Fontsize',FntSz+2)
hcb1.TickLabels{1} = '\leq0.05';

axis([minFA maxFA minFA maxFA minTR maxTR])
xticks([5 10 15 20]); xlabel('\alpha_2 (deg)')
yticks([5 10 15 20]); ylabel('\alpha_1 (deg)')
zticks([5 10 15 20]); zlabel('TR (ms)')

set(gca,'Fontsize',FntSz)


hcb1.Position(1) = hcb1.Position(1)+0.01;
hcb1.Position(2) = hcb1.Position(2)+0.15;
hcb1.Position(4) = hcb1.Position(4)-0.3;

hsp1.Position(1) = hsp1.Position(1)-0.01;
hsp1.Position(2) = hsp1.Position(2)+0.04;
hsp1.Position(3) = hsp1.Position(3)+0.02;
hsp1.Position(4) = hsp1.Position(4)-0.05;

set(gca,'SortMethod','ChildOrder');

hspFA = subplot(1,3,3);
text(-0.25,1.07,'c','Units','normalized','Fontsize',FntSz+6,'FontWeight','bold')

nFA = 100;
FAlist = linspace(0, pi/9, nFA);
sigSPGR = zeros(nFA, 1); 
for ii=1:nFA
    sigSPGR(ii) = signal_SPGR([FAlist(ii) uopt(3)],[T1 T2 M0]);
end
sigA1 = signal_SPGR([uopt(1) uopt(3)],[T1 T2 M0]);
sigA2 = signal_SPGR([uopt(2) uopt(3)],[T1 T2 M0]);

hold on;
plot(rad2deg(FAlist), sigSPGR, 'k', 'Linewidth', 2)
plot(rad2deg(uopt(1)), sigA1, 'o', 'Linewidth', 2, 'MarkerEdgeColor', c, 'MarkerFaceColor', c, 'MarkerSize', 10)
plot(rad2deg(uopt(2)), sigA2, 'o', 'Linewidth', 2, 'MarkerEdgeColor', c, 'MarkerFaceColor', c, 'MarkerSize', 10)

grid minor; box on;
yticks(0:0.02:max_y); ylabel('|M_{xy}|')
xticks(0:5:20); xlabel('Flip angle (deg)')
set(gca,'Fontsize',FntSz)

text(rad2deg(uopt(1))-2, sigA2+0.0075, [char(945),' = ',num2str(rad2deg(uopt(1)),'%.1f'),'{}^\circ'],...
    'Fontsize',FntSz,'Color',c)
text(rad2deg(uopt(2))+0.8, sigA1, [char(945),' = ',num2str(rad2deg(uopt(2)),'%.1f'),'{}^\circ'],...
    'Fontsize',FntSz,'Color',c)

text(10, 0.015,['TR = ',num2str(uopt(3),'%.1f'),'ms'],'Fontsize',FntSz,'HorizontalAlignment','center',...
    'Color',c)

hspFA.Position(1) = hspFA.Position(1)+0.05;
hspFA.Position(2) = hspFA.Position(2)+0.06;
hspFA.Position(4) = hspFA.Position(4)-0.07;


% print('figure1','-dtiff','-r600')



