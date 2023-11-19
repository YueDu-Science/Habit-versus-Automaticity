%% Figure 2F and 2H
clear all;
clc;
close all;
% load data
addpath(genpath('D:\Project\Publication_Data_Code\Habit-versus-Automaticity\analysis'))
load Data_for_Analysis;

zhongguo_hong = [162, 39, 41]/255;
jincao_zi = [148,153,187]/255;
jingui = [238,159,67]/255;
mks = 8; lw = 1;

% plot 
phit_remap_figure = figure('name','phit');

set(gca,'TickDir','out');
set(gca,'fontsize',10,'FontWeight','normal')
hAx=gca;                    % create an axes
hAx.LineWidth=0.3;
xplot = 0.001:0.001:1.2;
Xplot = xplot*1000;
axis([0 max(xplot) 0 1]);
set(gca,'YTick',[0:0.2:1],'YTicklabel',[0:0.2:1]);
set(gca,'YaxisLocation','left');
set(gca,'XTick',[0:0.3:1.2],'XTicklabel',[0:300:1200]);
set(gcf,'color','w');
hold on

swap_online = DATA_NA.Online(1).Phit;
% 5 bad datasets; Figure S1
ind = find(DATA_NA.Online(1).Subject ~= 118 & DATA_NA.Online(1).Subject ~= 123 ...
    & DATA_NA.Online(1).Subject ~= 128 & DATA_NA.Online(1).Subject ~= 130 & DATA_NA.Online(1).Subject ~= 139);
swap = cell2mat(swap_online(ind));

% get fitting to calculate t_min
para_online = DATA_NA.Online(1).Para;
para_online = para_online(ind);
para = cell2mat(para_online);
remap_mu = para(:,1);

% get fitting y_hat to calculate t_min
ycdf_online = DATA_NA.Online(1).Ycdf;
ycdf_online = ycdf_online(ind);
ycdf = cell2mat(ycdf_online);

% get each type of response to remapped stimuli
temp = reshape(swap',4,length(Xplot),numel(swap)/(length(Xplot)*4));
remapped = temp(2,:,:); 
unchange = temp(1,:,:); 
habitual = temp(4,:,:);
other_error = temp(3,:,:);
clear temp;
       
%shadedErrorBar(xplot,nanmean(unchange,3),seNaN(unchange),{'-','color',cols(1,:,2)})
shadedErrorBar(xplot,nanmean(remapped,3),seNaN(remapped),{'-','color',jincao_zi})
shadedErrorBar(xplot,nanmean(habitual,3),seNaN(habitual),{'-','color',zhongguo_hong})
shadedErrorBar(xplot,nanmean(other_error,3),seNaN(other_error),{'-','color',jingui})

%f1 = plot(xplot,nanmean(unchange,3),'-','color',cols(1,:,2),'Markersize',mks,'linewidth',lw);
f2 = plot(xplot,nanmean(remapped,3),'-','color',jincao_zi,'Markersize',mks,'linewidth',lw);
f3 = plot(xplot,nanmean(habitual,3),'-','color',zhongguo_hong,'Markersize',mks,'linewidth',lw);
f4 = plot(xplot,nanmean(other_error,3),'-','color',jingui,'Markersize',mks,'linewidth',lw);
  
plot([0 1.2], [0.125 0.125],'k--','linewidth',0.5);
text(0.5,0.2,'chance','FontSize',10,'FontWeight','Normal')   

legend([f2, f3, f4],{'Correct','Habitual','Other-Error'},'Location','North','NumColumns',1,...
     'fontsize',10,'textcolor','k');
    legend('boxoff');

posMat = get(gca,'Position');
posMat(3) = 0.6;
posMat(4) = 0.6;
set(gca,'Position',posMat);

% set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 10, 7],...
%  'PaperUnits', 'centimeters', 'PaperSize', [10, 7])
% cd('D:\Project\Habit_Formation\AVMA_Size\analysis\Figure');
% set(gcf,'renderer','Painters');
% print(phit_remap_figure,'Exp1_phit_remap', '-depsc','-r600');

%% Figure 2F inset
%%% normalize habitual curve (rt < tmin) to chance level
% to get rid of response bias which may contaminate habit measurement
% see text methods and Figure 2F

% get t_min for each participant
for s = 1:size(habitual,3)
    tmp = habitual(1,:,s);
    chance = para(s,3);a
    threshold = chance + 0.05;
    ind = find(ycdf(s,:) >= threshold, 1);
    tmin(s) = ind;
end

% based on t_min, align the baseline probability to the theoreical chance
% 0.125
for i = 1:size(habitual,3)
   a = nanmean(habitual(1,1:tmin(s),i)) - 0.125;
   habitual(1,:,i) = habitual(1,:,i) - a;
end


%% permutation test to see if the probability of habitual resposne is higher than chance (without bias)
%% calculate based on t_min
% adjust habitual and other-error curves based t_min
% because each participant had different speeds so habit may happen in
% different RT interval; aligning with t_min can allow us to identify habit
% in RT interval relative to t_min;

% use this to run permutation test after adjusted with tmin and plot the
% adjusted SAT; 

% habit_adjust_swap = [];
% for s = 1:size(habitual,3)
%     tmp = habitual(1,:,s);
%     chance = para(s,3);
%     threshold = chance + 0.05;
%     ind = find(ycdf(s,:) >= threshold, 1);
%     if ind < 300
%         tmp_adjust = [nan(1, 300-ind), tmp(1: ind), tmp(ind+1: 1200 - 300 + ind)];
%     else
%         tmp_adjust = [tmp(ind - 300 + 1: ind), tmp(ind+1: 1200), nan(1, ind-300)];
%     end
%     habit_adjust_swap(1,:,s) = tmp_adjust;
%     tmin(s) = ind;
% end
% 
% habitual = habit_adjust_swap;

% in our paper, we plotted the original curves without aligning to t_min;
habit = reshape(habitual, 1200, numel(habitual)/1200);
chance = repmat(0.125, 1200, size(habit,2));
other = reshape(other_error, 1200, numel(other_error)/1200);
clear temp;
x = 1:1200;
[C_global] = permutation_test_2tailed(habit(x,:),chance(x,:),1,1);
% find x where significant results were found
ind = find(nanmean((habit(x,:)-chance(x,:))') > C_global(1,:) | nanmean((habit(x,:)-chance(x,:))') < C_global(2,:));
% since we use a sliding window of 100ms
% we extend the ind by -50 and +50
IND_swap = ind;
IND_swap = [414-50,611+50]; % this is besed on the orignal data without adjusted with tmin;
                            % the average tmin is 372ms, after taking it
                            % out, it is generally matched with the
                            % IND_swa below

%IND_swap = [373-50,552+50]; % this is after adjusted based on tmin;
%adjusted by 50 ms because we used a 100ms smoothing window
%tmin was set to 300 for plot.

% plot Figure 2F inset
phit_remap_figure = figure('name','phit_adjust');
set(gca,'TickDir','out');
set(gca,'fontsize',10,'FontWeight','normal')
hAx=gca;                    % create an axes
hAx.LineWidth=0.3;
xplot = 0.001:0.001:1.2;
Xplot = xplot*1000;
axis([0 max(xplot) 0 0.5]);
set(gca,'YTick',[0:0.2:1],'YTicklabel',[0:0.2:1]);
set(gca,'YaxisLocation','left');
set(gca,'XTick',[0:0.3:1.2],'XTicklabel',[0:300:1200]);
%use this is plot after adjusted based on tmin
% set(gca,'XTick',[0:0.3:1.2],'XTicklabel',{'t_{min} - 300', 't_{min}'...
%     't_{min} + 300', 't_{min} + 600', 't_{min} + 900'});
%xtickangle(30);
set(gcf,'color','w');
hold on

shadedErrorBar(xplot,nanmean(habitual,3),seNaN(habitual),{'-','color',zhongguo_hong})
f3 = plot(xplot,nanmean(habitual,3),'-','color',zhongguo_hong,'Markersize',mks,'linewidth',lw);

% plot RT interval where permutation test is significant
f4 = plot(IND_swap/1000,ones(length(IND_swap),1) - 0.5,'-','color','r','Markersize',mks,'linewidth',lw+0.5);

plot([0 1.2], [0.125 0.125],'k--','linewidth',0.5);
text(0.5,0.2,'chance','FontSize',10,'FontWeight','Normal')   

posMat = get(gca,'Position');
posMat(3) = 0.6;
posMat(4) = 0.6;
set(gca,'Position',posMat);

% set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 10, 7],...
%  'PaperUnits', 'centimeters', 'PaperSize', [10, 7])
% cd('D:\Project\Habit_Formation\AVMA_Size\analysis\Figure');
% set(gcf,'renderer','Painters');
% print(phit_remap_figure,'Exp1_habit_index', '-depsc','-r600');


%% Figure 2H
swap_online_hi = DATA_NA.Online(1).habit_index_chance;
% 5 bad datasets
ind = find(DATA_NA.Online(1).Subject ~= 118 & DATA_NA.Online(1).Subject ~= 123 ...
    & DATA_NA.Online(1).Subject ~= 128 & DATA_NA.Online(1).Subject ~= 130 & DATA_NA.Online(1).Subject ~= 139);
swap_hi = cell2mat(swap_online_hi(ind));

swap_hi = swap_hi(:,2) - swap_hi(:,1); % adjust so baseline index was around 0 otherwise it inflates the magnitude of habit

mks = 4; lw = 1;
bi_xi = [52,94,100]/255;
brown = [139,94,60]/255;
grey = [0.7,0.7,0.7];

Diff_figure = figure('name','Diff');
set(gca,'TickDir','out');
set(gca,'fontsize',10,'FontWeight','normal')
hAx=gca;                    % create an axes
hAx.LineWidth=0.5;
ylabel('Habit index','FontSize',10, 'FontWeight','normal');
axis([0.5 1.5 -0.2 0.5]);
set(gca,'YTick',[-0.1:0.1:0.5],'YTicklabel',[-0.1:0.1:0.5]);
set(gca,'YaxisLocation','left');
set(gca,'XTick',[]);
set(gcf,'color','w');
hold on

f1 = bar(1,nanmean(swap_hi),0.3,'FaceColor','w','EdgeColor','r','Linewidth',lw);
alpha(f1,.3)

xx = 1+randn(100,1)*0.1;
for s = 1:size(swap_hi,1)
    if swap_hi(s) > 0 % if no habit
       g = plot(xx(s),swap_hi(s),'o','color','r','markerfacecolor','r','Markersize',mks,'linewidth',lw-0.4)
       g.Color(4) = 0.5;
    else
        plot(xx(s),swap_hi(s),'o','color','r','markerfacecolor','w','Markersize',mks,'linewidth',lw-0.4)
    end
end

h = [1];
hE = errorbar(h',nanmean(swap_hi),nanstd(swap_hi)/sqrt(numel(swap_hi)),...
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
% print(Diff_figure,'Exp1_phit_habit_index_raw', '-depsc','-r600');
