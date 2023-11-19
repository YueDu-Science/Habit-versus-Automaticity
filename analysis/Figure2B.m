%% Fig 2B
clear all;
clc;
close all;
% load data
addpath(genpath('D:\Project\Publication_Data_Code\Habit-versus-Automaticity\analysis'))
load Data_for_Analysis;

mks = 8; lw = 1; ms = 2;
zhu_dan = [176,82,76]/255;
bi_xi = [52,94,100]/255;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RT_online = DATA_NA.Online(1).RT;
RT = cell2mat(RT_online);
% 5 bad datasets; Figure S1
ind = find(DATA_NA.Online(1).Subject ~= 118 & DATA_NA.Online(1).Subject ~= 123 ...
    & DATA_NA.Online(1).Subject ~= 128 & DATA_NA.Online(1).Subject ~= 130 & DATA_NA.Online(1).Subject ~= 139);
RT = RT(ind,:);
rt4 = RT(:,1); rt8 = RT(:,2);

%%% plot RT 
RT_figure = figure('name','RT_MU');

set(gca,'TickDir','out');
set(gca,'fontsize',10,'FontWeight','normal')
hAx=gca;                    % create an axes
hAx.LineWidth=0.3; 

ylabel('RT (ms)','FontSize',10, 'FontWeight','normal');
%set(gca,'XTick',[4, 5], 'XTickLabel', {'Size 4', 'Size8'},'FontSize',12, 'FontWeight','normal');
set(gca,'XTick',[]);
axis([3 6 0.25 1])
ax = gca;
ax.XAxis.TickLength = [0,0]; 
set(gcf,'color','w');
hold on

f1 = bar(4,nanmean(rt4),0.7,'FaceColor',zhu_dan,'EdgeColor','None','Linewidth',lw);
f2 = bar(5,nanmean(rt8),0.7,'FaceColor',bi_xi,'EdgeColor','None','Linewidth',lw);

h = [4];
hE = errorbar(h',nanmean(rt4),nanstd(rt4)/sqrt(numel(rt4)),...
    'k');
set(hE(1),'LineWidth',0.5,'color','k')
hE.CapSize = 0;
% 
h = [5];
hE = errorbar(h',nanmean(rt8),nanstd(rt8)/sqrt(numel(rt8)),...
    'k');
set(hE(1),'LineWidth',0.5,'color','k')
hE.CapSize = 0;

for s = 1:size(rt4,1)
    plot([4.2 4.8],[rt4 rt8],'-k','color', qingshui_lan,...
        'markerfacecolor','w','Markersize',ms,'linewidth',0.3)
end

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
% print(RT_figure,'Exp1_rt_remap_bar', '-depsc','-r600');
% 
% [h,p,ci,stats] = ttest(rt4, rt8)
% cohend = nanmean(rt4*1000 - rt8*1000) / nanstd(rt4*1000 - rt8*1000);
% 
% diff_rt = nanmean(rt4*1000 - rt8*1000);
% std_rt = nanstd(rt4*1000 - rt8*1000)/sqrt(numel(rt4));
% 
% text(3.2,1,['\delta = ' num2str(diff_rt) '\pm' num2str(std_rt) ' ms'],'FontSize',10)


%% Confirm no sacrafice of accuracy 
AC_online = DATA_NA.Online(1).AC;
AC = cell2mat(AC_online);
ind = find(DATA_NA.Online(1).Subject ~= 118 & DATA_NA.Online(1).Subject ~= 123 ...
    & DATA_NA.Online(1).Subject ~= 128 & DATA_NA.Online(1).Subject ~= 130 & DATA_NA.Online(1).Subject ~= 139);
AC = AC(ind,:);

ac4 = AC(:,1); ac8 = AC(:,2);

[h,p,ci,stats] = ttest(ac4, ac8)
[p,h,stats] = ranksum(ac4, ac8)

diff_ac = nanmean(ac4*1000 - ac8*1000); % 4-element task had better accuracy
std_ac = nanstd(ac4*1000 - ac8*1000)/sqrt(numel(ac4));


