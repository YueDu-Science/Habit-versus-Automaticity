%% Fig 2I
clear all;
clc;
close all;
% load data
addpath(genpath('C:\Science\Publication_Data_Code\Habit-versus-Automaticity\analysis'))
load Data_for_Analysis;
% initialize parameters
grey = [150,150,150]/255;
brown = [244,164,96]/255;
mks = 8; lw = 1; 

% Habit index
swap_online_hi = DATA_NA.Online(1).habit_index;
% 5 bad datasets; Figure S1
ind = find(DATA_NA.Online(1).Subject ~= 118 & DATA_NA.Online(1).Subject ~= 123 ...
    & DATA_NA.Online(1).Subject ~= 128 & DATA_NA.Online(1).Subject ~= 130 & DATA_NA.Online(1).Subject ~= 139);
swap_hi = cell2mat(swap_online_hi(ind));
swap_hi = swap_hi(:,2) - swap_hi(:,1); % adjust so baseline index was around 0

% set size effect from SAT
para_online = DATA_NA.Online(1).SAT_Para;
para = cell2mat(para_online);
para = para(ind,:);
mu4 = para(:,1); mu8 = para(:,5);
mu = mu4 - mu8;

% subjects with ot without habit-determined by model comparison
ind_no = [19,22,23];
ind_h = [1:18,20,21];

% plot
SAT_MU_figure = figure('name','Diff');
set(gca,'TickDir','out');
set(gca,'fontsize',10,'FontWeight','normal')
hAx=gca;                    % create an axes
hAx.LineWidth=0.5;

%title('Withholding','FontSize',12, 'FontWeight','normal');
ylabel('\delta_{4-8}(ms)','FontSize',10, 'FontWeight','normal');
xlabel('Habit index','FontSize',10, 'FontWeight','normal');
axis([-0.12 0.55 -0.15 0.02]);
set(gca,'YTick',[-0.15:0.05:0],'YTicklabel',[-0.15:0.05:0]*1000);
set(gca, 'YDir','reverse')
set(gca,'YaxisLocation','left');
set(gca,'XTick',[-0.1:0.1:0.5],'XTicklabel',[-0.1:0.1:0.5]);
set(gcf,'color','w');
hold on

% data from habit people
plot(swap_hi(ind_h),mu(ind_h),'o','color',brown,'markerfacecolor',brown,'Markersize',mks,'linewidth',lw-0.4)

% correlation
[rho, p] = corr(swap_hi(ind_h), mu(ind_h))
text(0,-0.15,['\rho = ' num2str(rho) '; p = ' num2str(ceil(p*100)/100)],'FontSize', 10, 'color','b');
lsline

% data from no habit people
plot(swap_hi(ind_no),mu(ind_no),'o','markeredgecolor',brown,'markerfacecolor','w','Markersize',mks,'linewidth',lw-0.4)

posMat = get(gca,'Position');
posMat(4) = 0.5;
posMat(3) = 0.5;
set(gca,'Position',posMat);        
set(gca, 'Layer', 'top');
             
% set(gcf, 'Units', 'Centimeters', 'Position', [0, 0, 10, 7],...
%      'PaperUnits', 'Centimeters', 'PaperSize', [10, 7])
%  set(gcf,'renderer','Painters');
% cd('C:\Science\Project\Habit\AVMA_Size\analysis\Figure');
% set(gcf,'renderer','Painters');
% print(SAT_MU_figure,'Exp1_sat_vs_habit', '-depsc','-r600');

