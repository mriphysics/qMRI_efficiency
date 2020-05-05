%% Analysis of the under-sampling impact using dynamics factor dR
%David Leitao (david.leitao@kcl.ac.uk); 23-04-20

clearvars; close all; clc;

%% Load aliasing results from Monte Carlo simulations

% load random sampling results and save in structure random
load('random_sampling_SaR.mat')
random.SaR = abs(aliasing.avg./aliasing.std);
load('random_sampling_SNR30.mat')
random.SNR30 = abs(aliasing.avg./aliasing.std);
load('random_sampling_SNR50.mat')
random.SNR50 = abs(aliasing.avg./aliasing.std);
load('random_sampling_SNR80.mat')
random.SNR80 = abs(aliasing.avg./aliasing.std);
random.R = R;

% load spiral sampling results and save in structure spiral
load('spiral_sampling_SaR.mat')
spiral.SaR = abs(aliasing.avg./aliasing.std);
load('spiral_sampling_SNR30.mat')
spiral.SNR30 = abs(aliasing.avg./aliasing.std);
load('spiral_sampling_SNR50.mat')
spiral.SNR50 = abs(aliasing.avg./aliasing.std);
load('spiral_sampling_SNR80.mat')
spiral.SNR80 = abs(aliasing.avg./aliasing.std);
spiral.R = R;

ytrue = phantom(N); %gold standard signal

%% Calculate d-factor maps

% d-factor calculation using left side of equation 7 in the paper
random.dfactor_SNR30 = zeros(N, N, numel(random.R));
random.dfactor_SNR50 = zeros(N, N, numel(random.R));
random.dfactor_SNR80 = zeros(N, N, numel(random.R));
% d-factor calculation using equation 8 in the paper
random.dfactor_SNR30_equiv  = zeros(N, N, numel(random.R));
random.dfactor_SNR50_equiv  = zeros(N, N, numel(random.R));
random.dfactor_SNR80_equiv  = zeros(N, N, numel(random.R));
% calculate d-factor maps for random sampling
for cc=1:numel(random.R)
    SNR_R = random.SNR30(:,:,cc);
    sigma_noise = 1/30;
    SNRimg = ytrue / sigma_noise;
    random.dfactor_SNR30(:,:,cc) = ( 1/sqrt(random.R(cc)) ) * SNRimg ./ SNR_R;
    random.dfactor_SNR30_equiv(:,:,cc) = sqrt(1 + (SNRimg.^2) ./ (random.R(cc)*random.SaR(:,:,cc).^2) );

    SNR_R = random.SNR50(:,:,cc);
    sigma_noise = 1/50;
    SNRimg = ytrue / sigma_noise;
    random.dfactor_SNR50(:,:,cc) = ( 1/sqrt(random.R(cc)) ) * SNRimg ./ SNR_R;
    random.dfactor_SNR50_equiv(:,:,cc) = sqrt(1 + (SNRimg.^2) ./ (random.R(cc)*random.SaR(:,:,cc).^2) );
    
    SNR_R = random.SNR80(:,:,cc);
    sigma_noise = 1/80;
    SNRimg = ytrue / sigma_noise;
    random.dfactor_SNR80(:,:,cc) = ( 1/sqrt(random.R(cc)) ) * SNRimg ./ SNR_R;
    random.dfactor_SNR80_equiv(:,:,cc) = sqrt(1 + (SNRimg.^2) ./ (random.R(cc)*random.SaR(:,:,cc).^2) );
end


% d-factor calculation using left side of equation 7 in the paper
spiral.dfactor_SNR30 = zeros(N, N, numel(spiral.R));
spiral.dfactor_SNR50 = zeros(N, N, numel(spiral.R));
spiral.dfactor_SNR80 = zeros(N, N, numel(spiral.R));
% d-factor calculation using equation 8 in the paper
spiral.dfactor_SNR30_equiv  = zeros(N, N, numel(spiral.R));
spiral.dfactor_SNR50_equiv  = zeros(N, N, numel(spiral.R));
spiral.dfactor_SNR80_equiv  = zeros(N, N, numel(spiral.R));
% calculate d-factor maps for spiral sampling
for cc=1:numel(spiral.R)
    SNR_R = spiral.SNR30(:,:,cc);
    sigma_noise = 1/30;
    SNRimg = ytrue / sigma_noise;
    spiral.dfactor_SNR30(:,:,cc) = ( 1/sqrt(spiral.R(cc)) ) * SNRimg ./ SNR_R;
    spiral.dfactor_SNR30_equiv(:,:,cc) = sqrt(1 + (SNRimg.^2) ./ (spiral.R(cc)*spiral.SaR(:,:,cc).^2) );

    SNR_R = spiral.SNR50(:,:,cc);
    sigma_noise = 1/50;
    SNRimg = ytrue / sigma_noise;
    spiral.dfactor_SNR50(:,:,cc) = ( 1/sqrt(spiral.R(cc)) ) * SNRimg ./ SNR_R;
    spiral.dfactor_SNR50_equiv(:,:,cc) = sqrt(1 + (SNRimg.^2) ./ (spiral.R(cc)*spiral.SaR(:,:,cc).^2) );
    
    SNR_R = spiral.SNR80(:,:,cc);
    sigma_noise = 1/80;
    SNRimg = ytrue / sigma_noise;
    spiral.dfactor_SNR80(:,:,cc) = ( 1/sqrt(spiral.R(cc)) ) * SNRimg ./ SNR_R;
    spiral.dfactor_SNR80_equiv(:,:,cc) = sqrt(1 + (SNRimg.^2) ./ (spiral.R(cc)*spiral.SaR(:,:,cc).^2) );
end

%% Plot d-factor maps for random and spiral sampling
    
figure;
set(gcf,'Units','normalized','Outerposition',[0 0 1 1],'Color','w')

subplot(6,1,1)
cat_random_dfactor_SNR30 = [];
for rr=1:numel(random.R)
    cat_random_dfactor_SNR30 = cat(2,cat_random_dfactor_SNR30,random.dfactor_SNR30(:,:,rr));
end
imagesc(cat_random_dfactor_SNR30,[0 8]); hcb = colorbar;
xticks([]); yticks([])
text(-0.17,0.5,'SNR_{image}= 30','Units','normalized','Fontsize',20,'Fontweight','bold','Color','k')
for ii=1:numel(random.R)
    text(N/2+(ii-1)*N,-20,['R=',num2str(random.R(ii))],'Fontsize',18,'Fontweight','bold','HorizontalAlignment','center','Color','k')
end
colormap('inferno'); set(gca,'visible','off')
set(gca,'Fontsize',18)
hcb.Ticks = [0 4 8]; hcb.Color = 'k';
axis image
text(-0.17,1.5,'a','Fontsize',28,'Fontweight','bold','HorizontalAlignment','center','Color','k','Units','normalized')

hnd2 = subplot(6,1,2);
cat_dfactor_SNR50 = [];
for rr=1:numel(random.R)
    cat_dfactor_SNR50 = cat(2,cat_dfactor_SNR50,random.dfactor_SNR50(:,:,rr));
end
imagesc(cat_dfactor_SNR50,[0 14]); hcb = colorbar;
xticks([]); yticks([])
text(-0.17,0.5,'SNR_{image}= 50','Units','normalized','Fontsize',20,'Fontweight','bold','Color','k')
colormap('inferno'); hcb.Ticks = [0 7 14]; set(gca,'visible','off');  hcb.Color = 'k';
set(gca,'Fontsize',18)
axis image
hnd2.Position(2) = hnd2.Position(2)+0.01;

hnd3 = subplot(6,1,3);
cat_dfactor_SNR80 = [];
for rr=1:numel(random.R)
    cat_dfactor_SNR80 = cat(2,cat_dfactor_SNR80,random.dfactor_SNR80(:,:,rr));
end
imagesc(cat_dfactor_SNR80,[0 22]); hcb = colorbar;
xticks([]); yticks([])
text(-0.17,0.5,'SNR_{image}= 80','Units','normalized','Fontsize',20,'Fontweight','bold','Color','k')
colormap('inferno'); hcb.Ticks = [0 11 22]; set(gca,'visible','off');  hcb.Color = 'k';
set(gca,'Fontsize',18)
axis image
hnd3.Position(2) = hnd3.Position(2)+0.02;

hnd4 = subplot(6,1,4);
cat_random_dfactor_SNR30 = [];
for rr=1:numel(spiral.R)
    cat_random_dfactor_SNR30 = cat(2,cat_random_dfactor_SNR30,spiral.dfactor_SNR30(:,:,rr));
end
imagesc(cat_random_dfactor_SNR30,[0 8]); hcb = colorbar;
xticks([]); yticks([])
text(-0.17,0.5,'SNR_{image}= 30','Units','normalized','Fontsize',20,'Fontweight','bold','Color','k')
for ii=1:numel(spiral.R)
    text(N/2+(ii-1)*N,-20,['R=',num2str(spiral.R(ii))],'Fontsize',18,'Fontweight','bold','HorizontalAlignment','center','Color','k')
end
colormap('inferno')
set(gca,'Fontsize',18); set(gca,'visible','off');  hcb.Color = 'k';
axis image
hcb.Ticks = [0 4 8];
text(-0.17,1.5,'b','Fontsize',28,'Fontweight','bold','HorizontalAlignment','center','Color','k','Units','normalized')

hnd5 = subplot(6,1,5);
cat_dfactor_SNR50 = [];
for rr=1:numel(spiral.R)
    cat_dfactor_SNR50 = cat(2,cat_dfactor_SNR50,spiral.dfactor_SNR50(:,:,rr));
end
imagesc(cat_dfactor_SNR50,[0 12]); hcb = colorbar;
xticks([]); yticks([])
text(-0.17,0.5,'SNR_{image}= 50','Units','normalized','Fontsize',20,'Fontweight','bold','Color','k')
colormap('inferno'); hcb.Ticks = [0 6 12];
set(gca,'Fontsize',18); set(gca,'visible','off');  hcb.Color = 'k';
axis image

hnd6 = subplot(6,1,6);
cat_dfactor_SNR80 = [];
for rr=1:numel(spiral.R)
    cat_dfactor_SNR80 = cat(2,cat_dfactor_SNR80,spiral.dfactor_SNR80(:,:,rr));
end
imagesc(cat_dfactor_SNR80,[0 18]); hcb  = colorbar;
xticks([]); yticks([])
text(-0.17,0.5,'SNR_{image}= 80','Units','normalized','Fontsize',20,'Fontweight','bold','Color','k')
colormap('inferno'); hcb.Ticks = [0 9 18];
set(gca,'Fontsize',18); set(gca,'visible','off');  hcb.Color = 'k';
axis image

hnd6.Position(2) = hnd6.Position(2)-0.08;
hnd5.Position(2) = hnd5.Position(2)-0.08;
hnd4.Position(2) = hnd4.Position(2)-0.08;

% print('SI_figureS3','-dtiff','-r200')

%% Calculate average and standard deviation of d-factor maps

nonzero_mask = ytrue>0;
random.avg_dfactor = zeros(3, numel(random.R));
random.std_dfactor = zeros(3, numel(random.R));
for cc=1:numel(random.R)
    aux = random.dfactor_SNR30(:,:,cc);
    random.avg_dfactor(1,cc) = mean(aux(nonzero_mask(:)));
    random.std_dfactor(1,cc) = std(aux(nonzero_mask(:)));
    aux = random.dfactor_SNR50(:,:,cc);
    random.avg_dfactor(2,cc) = mean(aux(nonzero_mask(:)));
    random.std_dfactor(2,cc) = std(aux(nonzero_mask(:)));
    aux = random.dfactor_SNR80(:,:,cc);
    random.avg_dfactor(3,cc) = mean(aux(nonzero_mask(:)));
    random.std_dfactor(3,cc) = std(aux(nonzero_mask(:)));
end

spiral.avg_dfactor = zeros(3, numel(spiral.R));
spiral.std_dfactor = zeros(3, numel(spiral.R));
for cc=1:numel(spiral.R)
    aux = spiral.dfactor_SNR30(:,:,cc);
    spiral.avg_dfactor(1,cc) = mean(aux(nonzero_mask(:)));
    spiral.std_dfactor(1,cc) = std(aux(nonzero_mask(:)));
    aux = spiral.dfactor_SNR50(:,:,cc);
    spiral.avg_dfactor(2,cc) = mean(aux(nonzero_mask(:)));
    spiral.std_dfactor(2,cc) = std(aux(nonzero_mask(:)));
    aux = spiral.dfactor_SNR80(:,:,cc);
    spiral.avg_dfactor(3,cc) = mean(aux(nonzero_mask(:)));
    spiral.std_dfactor(3,cc) = std(aux(nonzero_mask(:)));
end

%% Plot average of d-factor maps

%For plotting
SNRlevels = [30 50 80];
Rvalues   = [2 4 12];

% hypothesis: std_aliasing \propto R^0.5
%       y = k * (x-1)^0.5   = k * (z)^0.5, with z = x - 1
% => LOGy = LOGk + 0.5*LOGz
%    LOGy - 0.5*LOGz = 1 * LOGb <=> LOGb = 1 \ [LOGy - 0.5*LOGz]

random.avg_aSR = zeros(numel(random.R), 1);
for rr=1:numel(random.R)
    aux = 1./random.SaR(:,:,rr);
    random.avg_aSR(rr) = mean(aux(nonzero_mask));
end
LOGy = log(random.avg_aSR(2:end));
LOGx = log(random.R(2:end)-1)';
coeff = [ones(numel(random.R)-1,1)] \ (LOGy-0.5*LOGx);
random.k = exp(coeff(1));

spiral.avg_aSR = zeros(numel(spiral.R), 1);
for rr=1:numel(spiral.R)
    aux = 1./spiral.SaR(:,:,rr);
    spiral.avg_aSR(rr) = mean(aux(nonzero_mask));
end
LOGy = log(spiral.avg_aSR(2:end));
LOGx = log(spiral.R(2:end)-1)';
coeff = [ones(numel(spiral.R)-1,1)] \ (LOGy-0.5*LOGx);
spiral.k = exp(coeff(1));


c1 = lines(6);
FntSz = 16;

figure;
set(gcf,'Units','normalized','Outerposition',[0.1 0 0.8 1],'Color','w')

hnd1 = subplot(2,3,1);
hold on;
plot(100:110, 100:110,'-','Color','k','Linewidth',2,'MarkerFaceColor','k')
plot(100:110, 100:110,'--','Color','k','Linewidth',2,'MarkerFaceColor','k')

plot(random.R, random.avg_dfactor(1,:),'-o','Color',c1(1,:),'Linewidth',2,'MarkerFaceColor',c1(1,:))
plot(random.R, random.avg_dfactor(2,:),'-o','Color',c1(2,:),'Linewidth',2,'MarkerFaceColor',c1(2,:))
plot(random.R, random.avg_dfactor(3,:),'-o','Color',c1(3,:),'Linewidth',2,'MarkerFaceColor',c1(3,:))

plot(spiral.R, spiral.avg_dfactor(1,:),'--o','Color',c1(1,:),'Linewidth',2,'MarkerFaceColor',c1(1,:))
plot(spiral.R, spiral.avg_dfactor(2,:),'--o','Color',c1(2,:),'Linewidth',2,'MarkerFaceColor',c1(2,:))
plot(spiral.R, spiral.avg_dfactor(3,:),'--o','Color',c1(3,:),'Linewidth',2,'MarkerFaceColor',c1(3,:))
box on; grid on;
set(gca,'Fontsize',FntSz)
axis([R(1) R(end) 0 20])
xticks(R)
ylabel('Average  d_R'); xlabel('Under-sampling factor R')
text(-0.2,1.05,'a','Units','normalized','Fontsize',24,'Fontweight','bold')
hlg = legend('Random','Spiral','location','northeast');

hnd1.Position(1) = 0.1;
hnd1.Position(2) = 0.58;
hnd1.Position(3) = 0.4;
hnd1.Position(4) = 0.37;

hlg.Position(1) = 0.28;
hlg.Position(2) = 0.872;


hnd2 = subplot(2,3,3);
hold on;
ii = 4;
plot(110,110,'-', 'Color','k','Linewidth',2,'MarkerFaceColor','k')
plot(110,110,'--', 'Color','k','Linewidth',2,'MarkerFaceColor','k')
for rr=Rvalues
    plot(SNRlevels, random.avg_dfactor(:,random.R==rr),'-o', 'Color',c1(ii,:),'Linewidth',2,'MarkerFaceColor',c1(ii,:))
    plot(SNRlevels, spiral.avg_dfactor(:,round(spiral.R)==rr),'--o','Color',c1(ii,:),'Linewidth',2,'MarkerFaceColor',c1(ii,:))
    ii = ii + 1;
end

box on; grid on;
set(gca,'Fontsize',FntSz)
axis([25 85 0 20])
xticks(SNRlevels)
ylabel('Average  d_R'); xlabel('SNR_{image}')
text(-0.2,1.05,'b','Units','normalized','Fontsize',24,'Fontweight','bold')
hlg2 = legend('Random','Spiral','location','northwest');

hlg2.Position(1) = 0.7;

hnd2.Position(1) = 0.61;
hnd2.Position(2) = 0.58;
hnd2.Position(3) = 0.35;
hnd2.Position(4) = 0.37;


hnd_X1 = subplot(2,3,5);
hnd_X1.Position(1) = 1.1;
hnd_X1.Position(2) = 1.1;
hnd_X1.Position(3) = 1.1;
hnd_X1.Position(4) = 1.1;
hold on; 
plot(100,100,'o','Color',c1(1,:),'Linewidth',3,'MarkerFaceColor',c1(1,:))
plot(100,100,'o','Color',c1(2,:),'Linewidth',3,'MarkerFaceColor',c1(2,:))
plot(100,100,'o','Color',c1(3,:),'Linewidth',3,'MarkerFaceColor',c1(3,:))
hlg_X1 = legend('SNR_{image}= 30','SNR_{image}= 50','SNR_{image}= 80');
hlg_X1.Position(1) = 0.39;
hlg_X1.Position(2) = 0.84;
set(gca,'Fontsize',FntSz)


hnd_X2 = subplot(2,3,6);
hnd_X2.Position(1) = 3.1;
hnd_X2.Position(2) = 3.1;
hnd_X2.Position(3) = 3.1;
hnd_X2.Position(4) = 3.1;
hold on; c = lines(3);
plot(100,100,'o','Color',c1(4,:),'Linewidth',3,'MarkerFaceColor',c1(4,:))
plot(100,100,'o','Color',c1(5,:),'Linewidth',3,'MarkerFaceColor',c1(5,:))
plot(100,100,'o','Color',c1(6,:),'Linewidth',3,'MarkerFaceColor',c1(6,:))
hlg_X2 = legend('R=2','R=4','R=12');
hlg_X2.Position(1) = 0.6275;
hlg_X2.Position(2) = 0.855;
set(gca,'Fontsize',FntSz)

hnd3 = subplot(2,3,4);
hold on;
plot(random.R(2:end),random.avg_aSR(2:end),'^','Color',c1(1,:),'Linewidth',2)
plot(spiral.R(2:end),spiral.avg_aSR(2:end),'^','Color',c1(2,:),'Linewidth',2)
plot(2:48,random.k*sqrt((2:48)-1),'-.','Color',c1(1,:),'Linewidth',2)
plot(2:48,spiral.k*sqrt((2:48)-1),'-.','Color',c1(2,:),'Linewidth',2)

box on; grid on;
set(gca,'Fontsize',FntSz)
axis([1 48 0 8])
xticks(R)
ylabel('alias-to-signal ratio'); xlabel('Under-sampling factor R')
text(-0.2,1.05,'c','Units','normalized','Fontsize',24,'Fontweight','bold')
hnd_leg = legend('(SaR_{image}^{random})^{-1}','(SaR_{image}^{spiral})^{-1}',...
    'Location','northwest');

hnd3.Position(1) = 0.1;
hnd3.Position(2) = 0.08;
hnd3.Position(3) = 0.4;
hnd3.Position(4) = 0.37;

hnd_X3 = subplot(2,3,6);
hnd_X3.Position(1) = 15;
hnd_X3.Position(2) = 15;
hnd_X3.Position(3) = 0.2;
hnd_X3.Position(4) = 0.2;
plot(100:110,100:110,'-.','Color','k','Linewidth',2)
hndLeg_X3 = legend('$k \times \sqrt{R-1}$', 'location','northwest');
hndLeg_X3.FontSize = 16;
hndLeg_X3.Interpreter = 'latex';

hndLeg_X3.Position(1) = 0.235;
hndLeg_X3.Position(2) = 0.405;


% print('figure6','-dtiff','-r600')
