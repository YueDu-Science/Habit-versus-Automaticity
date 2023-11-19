%% Figure 2E
clear all;
clc;
close all;
% load data
addpath(genpath('D:\Project\Publication_Data_Code\Habit-versus-Automaticity\analysis'))
load Data_for_Analysis;

ms = 2; lw = 1;
yin_hong = [176,82,76]/255;
grey = [150,150,150]/255;
tai_lan = [23,54,97]/255;

% plot RT for correct responses to revised stimuli
CR_RT_figure = figure('name','CR_RT');

set(gca,'TickDir','out');
set(gca,'fontsize',10,'FontWeight','normal')
hAx=gca;                    % create an axes
hAx.LineWidth=0.3; 
ylabel('RT: Correct responses to remapped stimuli (ms)','FontSize',10, 'FontWeight','normal');
set(gca,'XTick',[]);
ax = gca;
ax.XAxis.TickLength = [0,0]; 
axis([0 3 0 1.8])
set(gca,'YTick',[0:0.4:1.6],'YTicklabel',[0:400:1600]);
set(gcf,'color','w');
hold on

swap_online = cell2mat(DATA_NA.Online(1).CR.rt_revise);
% 5 bad dataset; see text and Fig S1
ind = find(DATA_NA.Online(1).Subject ~= 118 & DATA_NA.Online(1).Subject ~= 123 ...
    & DATA_NA.Online(1).Subject ~= 128 & DATA_NA.Online(1).Subject ~= 130 & DATA_NA.Online(1).Subject ~= 139);
swap = swap_online(ind,:);

f1 = bar(1,nanmean(swap(:,1)),0.7,'FaceColor',grey,'EdgeColor','None','Linewidth',lw);
f2 = bar(2,nanmean(swap(:,2)),0.7,'FaceColor',tai_lan,'EdgeColor','None','Linewidth',lw);

xx1 = 1+randn(1,1)*0.05;
xx2 = 3+randn(1,1)*0.05;
for s = 1:size(swap,1)
    plot([1.1 1.9],[swap(s,1) swap(s,2)],'-k','color', [0.5,0.5,0.5],'markerfacecolor','w','Markersize',ms,'linewidth',.3)
end

h = [1];
hE = errorbar(h',nanmean(swap(:,1)),nanstd(swap(:,1))/sqrt(numel(swap(:,1))),...
    'k');
set(hE(1),'LineWidth',0.5,'color','k')
hE.CapSize = 0;
% 
h = [2];
hE = errorbar(h',nanmean(swap(:,2)),nanstd(swap(:,2))/sqrt(numel(swap(:,2))),...
    'k');
set(hE(1),'LineWidth',0.5,'color','k')
hE.CapSize = 0;

posMat = get(gca,'Position');
posMat(4) = 0.6;
posMat(3) = 0.3;
set(gca,'Position',posMat);
set(gca, 'Layer', 'top');

% set(gcf, 'Units', 'Centimeters', 'Position', [0, 0, 10, 7],...
%      'PaperUnits', 'Centimeters', 'PaperSize', [10, 7])
%  set(gcf,'renderer','Painters');
% cd('D:\Project\Habit_Formation\AVMA_Size\analysis\Figure');
% set(gcf,'renderer','Painters');
% print(CR_RT_figure,'Exp1_CR_Remap_RT', '-depsc','-r600');
% 
% [h,p,ci,stats] = ttest(swap(:,2), swap(:,1))
% cohend = nanmean(swap(:,2)*1000 - swap(:,1)*1000) / nanstd(swap(:,2)*1000 - swap(:,1)*1000);
% 
% diff_rt = nanmean(swap(:,2) - swap(:,1));
% std_rt = nanstd(swap(:,2)*1000 - swap(:,1)*1000)/sqrt(numel(swap(:,2)));
