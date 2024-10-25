%% load raw data
% Exp 1
clear; clc; close all;
path = ('C:\Science\Publication_Data_Code\Habit-versus-Automaticity\data\Exp1'); % the data folder path
folder = dir(path);
cd(path);

Pav_Data = table;
csv = dir('*.*csv');

for i = 1:size(char(csv.name),1)
    T = readtable(csv(i).name);
    clear tmp;
    
    tmp = T(:,{'OS','browser',...
        'symb_map', 'Remap_Pair_1','Remap_Pair_2','participant', 'date',...
        'stim_val','proposed_choice','key_num','stim_type','key','block_type',...
        'pre_rt','set_size','remap','repeat_count','trial_Count'...
        'grp_swap','block_num','prep_time',...
        'actual_press','rt','actual_choice','correct'}...
        );
  
    Pav_Data = [Pav_Data; tmp];
end

cd('C:\Science\Publication_Data_Code\Habit-versus-Automaticity\analysis');
datafname = ['AVMA_Size_Letter_Online'];
save(datafname, 'Pav_Data');

%% Load Exp2 data
clear; clc; close all;

path = ('C:\Science\Publication_Data_Code\Habit-versus-Automaticity\data\Exp2'); % the data folder path
folder = dir(path);
cd(path);

Pav_Data = table;
csv = dir('*.*csv');

for i = 1:size(char(csv.name),1)
    T = readtable(csv(i).name);
    clear tmp;
    
    tmp = T(:,{'OS','browser',...
        'session','symb_map', 'Remap_Pair_1','Remap_Pair_2','participant', 'date',...
        'stim_val','proposed_choice','key_num','stim_type','key','block_type',...
        'pre_rt','set_size','remap','repeat_count','trial_Count'...
        'grp_swap','block_num','prep_time',...
        'actual_press','rt','actual_choice','correct'}...
        );
  
    Pav_Data = [Pav_Data; tmp];
end

cd('C:\Science\Publication_Data_Code\Habit-versus-Automaticity\analysis');
datafname = ['AVMA_Size_Symbol_Online'];
save(datafname, 'Pav_Data');
%% Since this data set has different variables names from what we collected in lab settings
%% change variables to be the same so that we can combine two datasets for further analyses

%% load Exp2 data - Min
clear; clc; close all;

path = ('C:\Science\Publication_Data_Code\Habit-versus-Automaticity\data\Min'); % the data folder path
folder = dir(path);
cd(path);

Pav_Data = table;
csv = dir('*.*csv');

for i = 1:size(char(csv.name),1)
    T = readtable(csv(i).name);
    clear tmp;
    
    tmp = T(:,{'OS','browser',...
        'session','symb_map', 'Remap_Pair_1','Remap_Pair_2','participant', 'date',...
        'stim_val','proposed_choice','key_num','stim_type','key','block_type',...
        'pre_rt','set_size','remap','repeat_count','trial_Count'...
        'grp_swap','block_num','prep_time',...
        'actual_press','rt','actual_choice','correct'}...
        );
  
    Pav_Data = [Pav_Data; tmp];
end

cd('C:\Science\Publication_Data_Code\Habit-versus-Automaticity\analysis');
datafname = ['AVMA_Size_Symbol_Min_Online'];
save(datafname, 'Pav_Data');
%% Since this data set has different variables names from what we collected in lab settings
%% change variables to be the same so that we can combine two datasets for further analyses