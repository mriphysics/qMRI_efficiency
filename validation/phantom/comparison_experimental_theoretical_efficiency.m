%% Comparison of experimental and theoretical efficiency for MRF and conv.
%David Leitao (david.leitao@kcl.ac.uk); 22-08-19

clearvars; close all; clc;

%% Load theoretical/experimental efficiencies for CM and MRF

load('steadystate_theoretical_experimental_efficiency.mat')
steadystate_theo_eff = theoretical_eff;
steadystate_exp_eff  = experimental_eff;
%number of signals associated to each combination:
steadystate_nsignals  = zeros(size(combinations,1), 1); 
for ii=1:size(combinations,1)
    steadystate_nsignals(ii) = numel([combinations{ii,1}, combinations{ii,2}]);
end

load('transient_theoretical_experimental_efficiency')
transient_theo_eff = mean(all_theoretical_eff,3);
transient_exp_eff  = mean(all_experimental_eff,3);
%number of signals associated to each combination:
transient_nsignals  = sum(w,2);


%% Plot experimental and theoretical efficiencies for both experiments

max_x = 0.3;
max_y = 0.3;

figure; 
set(gcf,'Units','Normalized','Outerposition',[0.1 0.1 0.8 0.9],'Color','w')

%Plot data for conventional methods
c = viridis(numel(unique(steadystate_nsignals)));
hnd1 = subplot(2,2,1);
aux_x = [0 max_x]; aux_y = [0 max_y];
hold on; plot(linspace(aux_x(1),aux_x(2),100),linspace(aux_x(1),aux_x(2),100),'-k','Linewidth',2)
for nn=1:size(steadystate_theo_eff,1)
    plot(steadystate_theo_eff(nn,1),  steadystate_exp_eff(nn,1), 's', ...,
        'Color', c(steadystate_nsignals(nn)+1-min(steadystate_nsignals),:), ...
        'MarkerFaceColor', c(steadystate_nsignals(nn)+1-min(steadystate_nsignals),:))
end
ylabel('Experimental efficiency (s^{-1/2})')
xlabel('Theoretical efficiency (s^{-1/2})')
aux_title = title('T_1');
axis([aux_x(:); aux_y(:)])
xticks([0 0.1 0.2 0.3]); yticks([0 0.1 0.2 0.3])
xticklabels({'0', '0.1', '0.2', '0.3' }); yticklabels({'0', '0.1', '0.2', '0.3' })
set(gca,'Fontsize',14)
grid minor; axis square; box on;
drawnow 
hnd1_lines = findobj('Type','line'); legend([hnd1_lines(end)],{'y=x'},'Location','northwest')
aux_text = text(-0.55, 0.5, 'DESPOT/JSR', 'Units', 'normalized', 'Fontweight', 'bold', ...
    'Fontsize', aux_title.FontSize, 'Rotation', 90, 'HorizontalAlignment', 'center');
text(-0.4, 1.1, 'a', 'Units', 'normalized', 'Fontweight', 'bold', ...
    'Fontsize', aux_title.FontSize+4, 'HorizontalAlignment', 'center');


hnd2 = subplot(2,2,2);
aux_x = [0 max_x]; aux_y = [0 max_y];
%%% Add fitted and y=x lines
hold on; plot(linspace(aux_x(1),aux_x(2),100),linspace(aux_x(1),aux_x(2),100),'-k','Linewidth',2)
for nn=1:size(steadystate_theo_eff,1)
    plot(steadystate_theo_eff(nn,2),  steadystate_exp_eff(nn,2), 's', ...,
        'Color', c(steadystate_nsignals(nn)+1-min(steadystate_nsignals),:), ...
        'MarkerFaceColor', c(steadystate_nsignals(nn)+1-min(steadystate_nsignals),:))
end
ylabel('Experimental efficiency (s^{-1/2})')
xlabel('Theoretical efficiency (s^{-1/2})')
title('T_2')
axis([aux_x(:); aux_y(:)])
xticks([0 0.1 0.2 0.3]); yticks([0 0.1 0.2 0.3])
xticklabels({'0', '0.1', '0.2', '0.3' }); yticklabels({'0', '0.1', '0.2', '0.3' })
set(gca,'Fontsize',14)
grid minor; axis square; box on;
drawnow
hnd2_lines = findobj('Type','line'); legend([hnd2_lines(end-numel(hnd1_lines))],{'y=x'},'Location','northwest')
colormap('viridis');
aux_cb = colorbar;
aux_cb.Position(1) = aux_cb.Position(1) + 0.01;
aux_cb.TickLabels = {num2str(unique(steadystate_nsignals))};
aux_cb.FontSize = 14;
aux_cb_title = title(aux_cb, '#steady-states', 'Fontsize', 14);
aux_cb_title.Position(2) = 1.05*aux_cb_title.Position(2);
hnd2.Position(1) = hnd2.Position(1) - 0.05;
text(-0.4, 1.1, 'b', 'Units', 'normalized', 'Fontweight', 'bold', ...
    'Fontsize', aux_title.FontSize+4, 'HorizontalAlignment', 'center');


%Plot data for MRF methods
unq_MRF_sig = unique(transient_nsignals);
MRF_idx_sig = zeros(size(transient_nsignals)); 
for kk=1:numel(transient_nsignals)
    MRF_idx_sig(kk) = find(unq_MRF_sig == transient_nsignals(kk));
end
c = viridis(numel(unique(transient_nsignals)));

hnd3 = subplot(2,2,3);

aux_x = [0 max_x]; aux_y = [0 max_y];
%%% Add y=x line
hold on;
plot(linspace(aux_x(1),aux_x(2),100),linspace(aux_x(1),aux_x(2),100),'-k','Linewidth',2)

for nn=1:size(transient_theo_eff,1)
    plot(transient_theo_eff(nn,1),  transient_exp_eff(nn,1), 's', ...,
        'Color', c(MRF_idx_sig(nn),:), 'MarkerFaceColor', c(MRF_idx_sig(nn),:))
end

ylabel('Experimental efficiency (s^{-1/2})')
xlabel('Theoretical efficiency (s^{-1/2})')
aux_title = title('T_1');
axis([aux_x(:); aux_y(:)])
xticks([0 0.1 0.2 0.3]); yticks([0 0.1 0.2 0.3])
xticklabels({'0', '0.1', '0.2', '0.3' }); yticklabels({'0', '0.1', '0.2', '0.3' })
set(gca,'Fontsize',14)
grid minor; axis square; box on;
drawnow
hnd3_lines = findobj('Type','line');
legend([hnd3_lines(end-numel(hnd2_lines))],{'y=x'},'Location','northwest')

aux_text = text(-0.55, 0.5, 'non-DE spoiled MRF', 'Units', 'normalized', 'Fontweight', 'bold', ...
    'Fontsize', aux_title.FontSize, 'Rotation', 90, 'HorizontalAlignment', 'center');

text(-0.4, 1.1, 'c', 'Units', 'normalized', 'Fontweight', 'bold', ...
    'Fontsize', aux_title.FontSize+4, 'HorizontalAlignment', 'center');


hnd4 = subplot(2,2,4);

aux_x = [0 max_x]; aux_y = [0 max_y];
%%% Add fitted and y=x lines
hold on;
plot(linspace(aux_x(1),aux_x(2),100),linspace(aux_x(1),aux_x(2),100),'-k','Linewidth',2)

for nn=1:size(transient_theo_eff,1)
    plot(transient_theo_eff(nn,2),  transient_exp_eff(nn,2), 's', ...,
        'Color', c(MRF_idx_sig(nn),:), 'MarkerFaceColor', c(MRF_idx_sig(nn),:))
end

ylabel('Experimental efficiency (s^{-1/2})')
xlabel('Theoretical efficiency (s^{-1/2})')
title('T_2')
axis([aux_x(:); aux_y(:)])
xticks([0 0.1 0.2 0.3]); yticks([0 0.1 0.2 0.3])
xticklabels({'0', '0.1', '0.2', '0.3' }); yticklabels({'0', '0.1', '0.2', '0.3' })
set(gca,'Fontsize',14)
grid minor; axis square; box on;
drawnow
hnd4_lines = findobj('Type','line');
legend([hnd4_lines(end-numel(hnd3_lines))],{'y=x'},'Location','northwest')

colormap('viridis');
aux_cb = colorbar;
aux_cb.Position(1) = aux_cb.Position(1) + 0.01;
aux_cb.Ticks = linspace(0, 1, 8);
aux_cb.TickLabels = {['   '; ...
                      num2str(unq_MRF_sig(2)); ...
                      '   '; ...
                      num2str(unq_MRF_sig(4)); ...
                      '   '; ...
                      num2str(unq_MRF_sig(6)); ...
                      '   '; ...
                      num2str(unq_MRF_sig(8))]};
aux_cb.FontSize = 14;
aux_cb_title = title(aux_cb, '#time-points', 'Fontsize', 14);
aux_cb_title.Position(2) = 1.05*aux_cb_title.Position(2);

hnd4.Position(1) = hnd4.Position(1) - 0.05;

text(-0.4, 1.1, 'd', 'Units', 'normalized', 'Fontweight', 'bold', ...
    'Fontsize', aux_title.FontSize+4, 'HorizontalAlignment', 'center');


% print('figure4','-dtiff','-r600')
