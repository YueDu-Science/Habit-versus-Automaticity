%% Figure 2C and 2G
clear all;
clc;
close all;
% load data
addpath(genpath('D:\Project\Publication_Data_Code\Habit-versus-Automaticity\analysis'))
load Data_for_Analysis;

xplot = 1:1:1200;
mks = 8; lw = 1;

% color
zhu_dan = [176,82,76]/255;
bi_xi = [52,94,100]/255;
napoli_huang = [224,200,159]/255;
yin_hong = [176,82,76]/255;
tai_lan = [23,54,97]/255;

tiexiu_hong = [104, 20, 20]/255;
putaoyuan_lv = [79,89,62]/255;
qingshui_lan = [106,196,204]/255;

zhongguo_hong = [162, 39, 41]/255;
jincao_zi = [148,153,187]/255;
jingui = [238,159,67]/255;

% plot 
SAT_remap_figure = figure('name','SAT');
set(gca,'TickDir','out');
set(gca,'fontsize',10,'FontWeight','normal')
hAx=gca;                    % create an axes
hAx.LineWidth=0.3; 
xplot = 0.001:0.001:1.2;
Xplot = xplot*1000;
axis([0 max(xplot) -0.1 1.1]);
set(gca,'YTick',[0:0.2:1],'YTicklabel',[0:0.2:1]);
set(gca,'YaxisLocation','left');
set(gca,'XTick',[0:0.3:1.2],'XTicklabel',[0:300:1200]);
ylabel('Probability Correct','FontSize',10, 'FontWeight','normal');
xlabel('Allowed RT (ms)','FontSize',10, 'FontWeight','normal')
set(gcf,'color','w');
hold on
% sat data
swap_online = DATA_NA.Online(1).SAT;
% 5 bad datasets; figure S1
ind = find(DATA_NA.Online(1).Subject ~= 118 & DATA_NA.Online(1).Subject ~= 123 ...
    & DATA_NA.Online(1).Subject ~= 128 & DATA_NA.Online(1).Subject ~= 130 & DATA_NA.Online(1).Subject ~= 139);
swap = cell2mat(swap_online(ind));

% get fitting parameter \mu
para_online = DATA_NA.Online(1).SAT_Para;
para = cell2mat(para_online);
para = para(ind,:);
mu4 = para(:,1); mu8 = para(:,5);

% reorganize sat data for plot
temp = reshape(swap',2,length(Xplot),numel(swap)/(length(Xplot)*2));
sat4 = temp(1,:,:); 
sat8 = temp(2,:,:); 

% two different chance level
plot([0 max(Xplot)], [0.25 0.25],'k--','linewidth',0.5);
plot([0 max(Xplot)], [0.125 0.125],'k--','linewidth',0.5);
text(0.9,0.2,'chance','FontSize',12,'FontWeight','Normal')

% run this part to normalize SAT
% inset of Figure 2C
% % scaling so that all SAT ranges from 0 to 1;
% for i = 1:size(sat4,3)
%    sat4(1,:,i) = (sat4(1,:,i) - para(i,3))/(para(i,4) - para(i,3)); 
% end
% 
% for i = 1:size(sat8,3)
%    sat8(1,:,i) = (sat8(1,:,i) - para(i,7))/(para(i,8) - para(i,7)); 
% end
% 
% plot([0 max(Xplot)], [0 0],'k--','linewidth',0.5);

shadedErrorBar(xplot,nanmean(sat4,3),seNaN(sat4),{'-','color',zhu_dan})
shadedErrorBar(xplot,nanmean(sat8,3),seNaN(sat8),{'-','color',bi_xi})

f1 = plot(xplot,nanmean(sat4,3),'-','color',zhu_dan,'Markersize',mks,'linewidth',lw);
f2 = plot(xplot,nanmean(sat8,3),'-','color',bi_xi,'Markersize',mks,'linewidth',lw);

legend([f1, f2],{'4-element','8-element'},'Location','East','NumColumns',1,...
     'fontsize',10,'textcolor','k', 'FontWeight','normal');
             legend('boxoff');

posMat = get(gca,'Position');
posMat(3) = 0.6;
posMat(4) = 0.6;
set(gca,'Position',posMat);

% set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 10, 7],...
%      'PaperUnits', 'centimeters', 'PaperSize', [10, 7])
% cd('D:\Project\Habit_Formation\AVMA_Size\analysis\Figure');
% set(gcf,'renderer','Painters');
% print(SAT_remap_figure,'Exp1_sat_remap', '-depsc','-r600');

%% Figure 2C inset
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ms = 2; lw = 1;
%plot 
SAT_MU_figure = figure('name','SAT_MU');
set(gca,'TickDir','out');
set(gca,'fontsize',10,'FontWeight','normal')
hAx=gca;                    % create an axes
hAx.LineWidth=0.3; 
ylabel('Center (ms)','FontSize',10, 'FontWeight','normal');
set(gca,'XTick',[]);
axis([3 6 0.1 0.65])
set(gca,'YTick',[0:0.2:0.8],'YTicklabel',[0:200:800]);
ax = gca;
ax.XAxis.TickLength = [0,0]; 
set(gcf,'color','w');
hold on

f1 = bar(4,nanmean(mu4),0.7,'FaceColor',zhu_dan,'EdgeColor','None','Linewidth',lw);
f2 = bar(5,nanmean(mu8),0.7,'FaceColor',bi_xi,'EdgeColor','None','Linewidth',lw);

h = [4];
hE = errorbar(h',nanmean(mu4),nanstd(mu4)/sqrt(numel(mu4)),...
    'k');
set(hE(1),'LineWidth',0.5,'color','k')
hE.CapSize = 0;
% 
h = [5];
hE = errorbar(h',nanmean(mu8),nanstd(mu8)/sqrt(numel(mu8)),...
    'k');
set(hE(1),'LineWidth',0.5,'color','k')
hE.CapSize = 0;

for s = 1:size(mu4,1)
    plot([4.2 4.8],[mu4 mu8],'-k','color', qingshui_lan,...
        'markerfacecolor','w','Markersize',mks,'linewidth',.3)
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
% print(SAT_MU_figure,'Exp1_sat_remap_bar', '-depsc','-r600');
% 
% [h,p,ci,stats] = ttest(mu4, mu8)
% cohend = nanmean(mu4*1000 - mu8*1000) / nanstd(mu4*1000 - mu8*1000);
% 
% diff_mu = nanmean(mu4*1000 - mu8*1000);
% std_mu = nanstd(mu4*1000 - mu8*1000)/sqrt(numel(mu4));
% 
% 
% text(3.2,0.65,['\delta = ' num2str(diff_mu) '\pm' num2str(std_mu) ' ms'],'FontSize',10)

%% Figure 2G 
% plot individual delta
mks = 4; lw = 1;
swap = mu4 - mu8;

Diff_figure = figure('name','Diff');
set(gca,'TickDir','out');
set(gca,'fontsize',10,'FontWeight','normal')
hAx=gca;                    % create an axes
hAx.LineWidth=0.5;
ylabel('\delta_{4-8} (ms)','FontSize',10, 'FontWeight','normal');
axis([0.5 1.5 -0.16 0.015]);
set(gca,'YTick',[-0.3:0.05:0],'YTicklabel',[-0.3:0.05:0]*1000);
set(gca, 'YDir','reverse')
set(gca,'YaxisLocation','left');
set(gca,'XTick',[]);
set(gcf,'color','w');
hold on

f1 = bar(1,nanmean(swap),0.3,'FaceColor','w','EdgeColor',qingshui_lan,'Linewidth',lw);
alpha(f1,.3)

xx = 1+randn(100,1)*0.1;
for s = 1:size(swap,1)
    if s == 19 || s == 22 || s == 23 % participants who did not show habit; Figure 2H
        plot(xx(s),swap(s),'o','color',qingshui_lan,'markerfacecolor','w','Markersize',mks,'linewidth',lw-0.4)
    else
        g = plot(xx(s),swap(s),'o','color',qingshui_lan,'markerfacecolor',qingshui_lan,'Markersize',mks,'linewidth',lw-0.4)
        g.Color(4) = 0.5;
    end
end

h = [1];
hE = errorbar(h',nanmean(swap),nanstd(swap)/sqrt(numel(swap)),...
    'k','MarkerSize',3);
set(hE(1),'LineWidth',0.5,'color','k')
hE.CapSize = 0;

posMat = get(gca,'Position');
posMat(3) = 0.2;
posMat(4) = 0.8;
set(gca,'Position',posMat);

% set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 10, 7],...
%  'PaperUnits', 'centimeters', 'PaperSize', [10, 7])
% cd('D:\Project\Habit_Formation\AVMA_Size\analysis\Figure');
% set(gcf,'renderer','Painters');
% print(Diff_figure,'Exp1_sat_remap_bar_delta', '-depsc','-r600');
