%% Plot optimal results for the simple MRF experiment
%David Leitao (david.leitao@kcl.ac.uk); 07-01-21

clearvars; close all; clc;

%% Load optimal acquisition settings

load('opt_DE_MRF_5meas.mat')
N = 5;

%% Plot results for simple MRF experiments

figure;
set(gcf,'Units','normalized','Outerposition',[0.3 0.25 0.4 0.5],'Color','w')

LW = 2;     
dr     = 0.03; 

freq = 0.25;
fsinc =@(npts) 0.2 * cos(freq*linspace(-20,20,npts)) .* hanning(npts)';
c = lines(7);
cgray = [0.825 0.825 0.825];


%%% subplot (a)
hnd1 = subplot(2,1,1);
hnd1.Position = [0.05 0.6 0.85 0.35];
set(gca,'visible','off'); hold on; 
text(0,1,'(a)','Units','normalized','FontSize',18,'FontWeight','bold')
plot([-0.5 4.7],[0 0], '-k', 'Linewidth', LW)
xpos = [0.65 1.1 1.9 2.35 3.5]-0.5; dx = 0.15;
text(xpos(1)-0.4,0,'//', 'Fontsize', 14, 'HorizontalAlignment','center')
text(xpos(end)+0.4,0,'//', 'Fontsize', 14, 'HorizontalAlignment','center')

for ii=1:N
    sf = DE_joint.Acq(ii)/pi;
    if sf>0
        plot(linspace(xpos(ii)-dx,xpos(ii)+dx,1000),fsinc(1000)*sf,'Linewidth', LW, 'Color', 'k')
        ypos = max(fsinc(400))*sf;
    else
        ypos = 0;
        plot(xpos(ii)*ones(2,1),[-0.02 0.02],'-k','Linewidth',2)
    end

    if ii==1
        text(xpos(ii)+0.02, -6*dr, [num2str(sum(DE_joint.Acq(N+1:N+ii-1)),'%.0f'),'ms'], 'Fontsize', 12, 'HorizontalAlignment','center')
    elseif ii==3
        text(xpos(ii)-0.05, -6*dr, [num2str(sum(DE_joint.Acq(N+1:N+ii-1)),'%.0f'),'ms'], 'Fontsize', 12, 'HorizontalAlignment','center')
    elseif ii==4
        text(xpos(ii)+0.05, -6*dr, [num2str(sum(DE_joint.Acq(N+1:N+ii-1)),'%.0f'),'ms'], 'Fontsize', 12, 'HorizontalAlignment','center')
    else
        text(xpos(ii), -6*dr, [num2str(sum(DE_joint.Acq(N+1:N+ii-1)),'%.0f'),'ms'], 'Fontsize', 12, 'HorizontalAlignment','center')
    end
    text(xpos(ii), ypos+0.12, ['  ',num2str(rad2deg(DE_joint.Acq(ii)),'%.0f'),'\circ'], 'Fontsize', 12, 'HorizontalAlignment','center')
end

sf = DE_joint.Acq(1)/pi;
plot(linspace(xpos(end)+0.8-dx,xpos(end)+0.8+dx,400),fsinc(400)*sf,'Linewidth', LW, 'Color', 'k')
ypos = max(fsinc(400))*sf;
text(xpos(end)+0.8, ypos+0.12, [' ',num2str(rad2deg(DE_joint.Acq(1)),'%.0f'),'\circ'], 'Fontsize', 12, 'HorizontalAlignment','center')

sf = DE_joint.Acq(2)/pi;
plot(linspace(xpos(end)+0.8-dx+xpos(2)-xpos(1),xpos(end)+0.8+dx+xpos(2)-xpos(1),400),fsinc(400)*sf,'Linewidth', LW, 'Color', 'k')
ypos = max(fsinc(400))*sf;
text(xpos(end)+0.8+xpos(2)-xpos(1), ypos+0.12, [' ',num2str(rad2deg(DE_joint.Acq(2)),'%.0f'),'\circ'], 'Fontsize', 12, 'HorizontalAlignment','center')

fill([-0.5 xpos(1)-dx xpos(1)-dx -0.5],[0.5 0.5 -0.3 -0.3],'k','FaceAlpha',0.3,'EdgeAlpha',0,'LineStyle','none','EdgeColor','none')
fill([xpos(end)+0.8-dx 4.7 4.7 xpos(end)+0.8-dx],[0.5 0.5 -0.3 -0.3],'k','FaceAlpha',0.3,'EdgeAlpha',0,'LineStyle','none','EdgeColor','none')

plot([xpos(1)-dx xpos(end)+0.8-dx],[-0.3 -0.3],'k','LineWidth',LW)
plot(xpos(1)-dx+1.5*[dr 0 dr],1.2*[-dr 0 dr]-0.3, '-k', 'Linewidth', LW)
plot(xpos(end)+0.8-dx+1.5*[-dr 0 -dr],1.2*[-dr 0 dr]-0.3, '-k', 'Linewidth', LW)
text(mean([xpos(end)+0.8-dx xpos(1)-dx]), -0.375, [num2str(sum(DE_joint.Acq(N+1:N+6-1)),'%.0f'),'ms'], 'Fontsize', 12, 'HorizontalAlignment','center')


%%% subplot (b)
hnd2 = subplot(2,1,2);
hnd2.Position = [0.05 0.15 0.85 0.35];
set(gca,'visible','off'); hold on; 
text(0,1.1,'(b)','Units','normalized','FontSize',18,'FontWeight','bold')
plot([0 xpos(end)+0.8-dx],[0.5 0.5],'-','Linewidth',0.5,'Color',cgray); 
plot([0 0 xpos(end)+0.8-dx],[1 0 0], '-k', 'Linewidth', LW)
plot([0 xpos(end)+0.8-dx xpos(end)+0.8-dx],[1 1 0], '-k', 'Linewidth', 0.5)
text((xpos(end)+0.8-dx)/2,-14*dr,'t (ms)', 'Fontsize', 12, 'HorizontalAlignment','center')
plot([0 0.05],[0.5 0.5],'-k','Linewidth',2); 
text(-0.1,0,'0','FontSize',12,'HorizontalAlignment','right')
text(-0.1,0.5,'0.5','FontSize',12,'HorizontalAlignment','right')
text(-0.1,1,'1','FontSize',12,'HorizontalAlignment','right')

for ii=1:N    
    plot(xpos(ii)*ones(2,1),[-0.02 0.02],'-k','Linewidth',2)
    plot(xpos(ii),abs(DE_joint.Mxy(ii)),'o','Marker','o','MarkerFaceColor',c(1,:),'MarkerEdgeColor',c(1,:))
    plot(xpos(ii),781*abs(DE_joint.dMxydT1(ii)),'s','Marker','s','MarkerFaceColor',c(2,:),'MarkerEdgeColor',c(2,:))
    plot(xpos(ii),65*abs(DE_joint.dMxydT2(ii)),'^','Marker','^','MarkerFaceColor',c(3,:),'MarkerEdgeColor',c(3,:))
    text(xpos(ii), -6*dr,num2str(sum(DE_joint.Acq(N+1:N+ii-1))+2,'%.0f'), 'Fontsize', 12, 'HorizontalAlignment','center')
end
plot(-0.5,-0.3,'s','Color','w')
plot(4.7,1,'s','Color','w')

hlines = findobj(hnd2,'Type', 'line');
hlegend = legend(hlines(5:-1:3),'s','T_1\times\partials/\partialT_1','T_2\times\partials/\partialT_2',...
    'Orientation','Vertical','Location','best');
hlegend.FontSize = 11;
hlegend.Position = [0.768 0.363 0.06 0.12];


% print('figure1','-dtiff','-r600')
% print('figure1','-deps')
