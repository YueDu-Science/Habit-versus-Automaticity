%% load raw data

clear; clc; close all;
path = ('D:\Project\Publication_Data_Code\Habit-versus-Automaticity\data'); % the data folder path
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

cd('D:\Project\Publication_Data_Code\Habit-versus-Automaticity\analysis');
datafname = ['AVMA_Size_Letter_Online'];
save(datafname, 'Pav_Data');