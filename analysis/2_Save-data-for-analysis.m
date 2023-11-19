%% create a table file containing data that can be used for analysis and plot
clear all;
clc;
addpath(genpath('D:\Project\Publication_Data_Code\Habit-versus-Automaticity\analysis'))

%% load online data
load AVMA_Size_Letter_Online.mat;  % D
DATA = Pav_Data;
%head(DATA,3)  % check if data look right
max_rt = 1.3;
prac_blk = 1:10;
num_trg = 8;
xplot = 1:1:1200;
t_size = 100; %ms
group = 1; % only one group; within subject design

for g = 1:group
    tmp_sub = DATA.participant;
    sub_name = unique(tmp_sub);
    % deal with individual subject's data
    for s = 1:length(sub_name)
        data = [];
        data = DATA(DATA.participant == sub_name(s),:);
        
        Data_online(g).Subject(s) = sub_name(s);
        responsedata(g).Subject(s) = sub_name(s); % this dataset is used for modeling; not included in this paper
        % add a column to indicate handedness
        % these four subjects are left-handed
        if sub_name(s) == 127 || sub_name(s) == 138 || sub_name(s) == 139 || sub_name(s) == 140 
            data.handedness(:) = 'l';
        else
            data.handedness(:) = 'r';
        end
        
        % get remaped stimuli/responses
        stim_map = data.symb_map(~cellfun(@isempty,data.symb_map));
        stim_map = str2num(char(stim_map));
        pair_1 = data.Remap_Pair_1(~cellfun(@isempty,data.Remap_Pair_1));
        pair_2 = data.Remap_Pair_2(~cellfun(@isempty,data.Remap_Pair_2));

        remap_pair = [str2num(char(pair_1)); str2num(char(pair_2))];

%% Training session; free-RT        
        ind_blk_type = (all(char(data.block_type) == 'RT',2));
        ind_stim_type = (all(char(data.stim_type) == 'Symb',2));
        ind_pre_rt = (data.pre_rt == 1); % indicate not set-size effect test
        data_rt_symb_prac = data(ind_blk_type == 1 & ind_stim_type == 1 & ind_pre_rt,:);
        
        % remove any RT > 3000 ms;
        data_rt_symb_prac(data_rt_symb_prac.rt > 3,:) = [];

        for b = prac_blk
            rt = data_rt_symb_prac.rt(data_rt_symb_prac.block_num == b & data_rt_symb_prac.correct == 1);
            error = data_rt_symb_prac.correct(data_rt_symb_prac.block_num == b);
            mean_prep_time(b) = nanmean(rt);
            mean_correct(b) = nansum(error);
        end
        
        % save data
        Data_online(g).Prac.rt{s,1} = mean_prep_time;
        Data_online(g).Prac.ac{s,1} = mean_correct;
%% Learning sessions; responding as slowly as they need          
        index_blk_type = (all(char(data.block_type) == 'CR',2));
       
        index_no_change = (data.remap == 0); % original mapping
        index_new = (data.remap == 1); % revised mapping
        index_remap = (data.grp_swap == 1); % not useful as all participants are in the swap group in this study

        data_cr_symb_old = data(index_blk_type == 1 & index_no_change == 1,:);
        data_cr_symb_remap = data(index_blk_type == 1 & index_remap == 1 & index_new == 1,:);
        
        % remove RT > 8000 ms;
        data_cr_symb_old(data_cr_symb_old.rt > 8,:) = [];
        data_cr_symb_remap(data_cr_symb_remap.rt > 8,:) = [];
        
        
        % # trials they need to learn
        N_map(g,s) = max(data_cr_symb_old.trial_Count);
        N_remap(g,s) = max(data_cr_symb_remap.trial_Count);

        % RT for all trials
        rt_map(g,s) = nanmean(data_cr_symb_old.rt(data_cr_symb_old.repeat_count == 0));
        rt_remap(g,s) = nanmean(data_cr_symb_remap.rt(data_cr_symb_remap.repeat_count == 0));

        Data_online(g).CR.N{s,1} = N_map(g,s); % save data for plot group mean later
        Data_online(g).CR.N{s,2} = N_remap(g,s); 
        Data_online(g).CR.rt{s,1} = rt_map(g,s); 
        Data_online(g).CR.rt{s,2} = rt_remap(g,s);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%% calculate RT for different types of trials %%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % RT when they made correct errors to revised stimuli when learning

            if ~isempty(data_cr_symb_old)
                for trial = 1:size(data_cr_symb_old,1)
                    if data_cr_symb_old.proposed_choice(trial) == remap_pair(1,1)
                        data_cr_symb_old.is_remapped(trial) = 1;
                        data_cr_symb_old.remapped_from(trial) = remap_pair(1,2);
                    elseif data_cr_symb_old.proposed_choice(trial) == remap_pair(1,2)
                        data_cr_symb_old.is_remapped(trial) = 1;
                        data_cr_symb_old.remapped_from(trial) = remap_pair(1,1);
                    elseif data_cr_symb_old.proposed_choice(trial) == remap_pair(2,1)
                        data_cr_symb_old.is_remapped(trial) = 1;
                        data_cr_symb_old.remapped_from(trial) = remap_pair(2,2);
                    elseif data_cr_symb_old.proposed_choice(trial) == remap_pair(2,2)
                        data_cr_symb_old.is_remapped(trial) = 1;
                        data_cr_symb_old.remapped_from(trial) = remap_pair(2,1);
                    else
                        data_cr_symb_old.is_remapped(trial) = 0;
                        data_cr_symb_old.remapped_from(trial) = nan;
                    end
                end
            end
            cr_habit_remap = data_cr_symb_old(data_cr_symb_old.is_remapped == 1,:);
            cr_habit_remap_trials = cr_habit_remap(cr_habit_remap.correct == 1,:);
            Data_online(g).CR.rt_revise{s,1} = nanmean(cr_habit_remap_trials.rt);
            Data_online(g).CR.N_revise{s,1} = length(cr_habit_remap_trials.rt);

            if ~isempty(data_cr_symb_remap)
                for trial = 1:size(data_cr_symb_remap,1)
                    if data_cr_symb_remap.proposed_choice(trial) == remap_pair(1,1)
                        data_cr_symb_remap.is_remapped(trial) = 1;
                        data_cr_symb_remap.remapped_from(trial) = remap_pair(1,2);
                    elseif data_cr_symb_remap.proposed_choice(trial) == remap_pair(1,2)
                        data_cr_symb_remap.is_remapped(trial) = 1;
                        data_cr_symb_remap.remapped_from(trial) = remap_pair(1,1);
                    elseif data_cr_symb_remap.proposed_choice(trial) == remap_pair(2,1)
                        data_cr_symb_remap.is_remapped(trial) = 1;
                        data_cr_symb_remap.remapped_from(trial) = remap_pair(2,2);
                    elseif data_cr_symb_remap.proposed_choice(trial) == remap_pair(2,2)
                        data_cr_symb_remap.is_remapped(trial) = 1;
                        data_cr_symb_remap.remapped_from(trial) = remap_pair(2,1);
                    else
                        data_cr_symb_remap.is_remapped(trial) = 0;
                        data_cr_symb_remap.remapped_from(trial) = nan;
                    end
                end
            end
            cr_habit_remap = data_cr_symb_remap(data_cr_symb_remap.is_remapped == 1,:);
            cr_habit_remap_trials = cr_habit_remap(cr_habit_remap.actual_choice == cr_habit_remap.proposed_choice | cr_habit_remap.actual_choice + 4 == cr_habit_remap.proposed_choice,:);
            Data_online(g).CR.rt_revise{s,2} = nanmean(cr_habit_remap_trials.rt);
            Data_online(g).CR.N_revise{s,2} = length(cr_habit_remap_trials.rt);

            
            % RT when they made correct errors to non-revised stimuli when learning
            if ~isempty(data_cr_symb_old)
                for trial = 1:size(data_cr_symb_old,1)
                    if data_cr_symb_old.proposed_choice(trial) == remap_pair(1,1)
                        data_cr_symb_old.is_remapped(trial) = 1;
                        data_cr_symb_old.remapped_from(trial) = remap_pair(1,2);
                    elseif data_cr_symb_old.proposed_choice(trial) == remap_pair(1,2)
                        data_cr_symb_old.is_remapped(trial) = 1;
                        data_cr_symb_old.remapped_from(trial) = remap_pair(1,1);
                    elseif data_cr_symb_old.proposed_choice(trial) == remap_pair(2,1)
                        data_cr_symb_old.is_remapped(trial) = 1;
                        data_cr_symb_old.remapped_from(trial) = remap_pair(2,2);
                    elseif data_cr_symb_old.proposed_choice(trial) == remap_pair(2,2)
                        data_cr_symb_old.is_remapped(trial) = 1;
                        data_cr_symb_old.remapped_from(trial) = remap_pair(2,1);
                    else
                        data_cr_symb_old.is_remapped(trial) = 0;
                        data_cr_symb_old.remapped_from(trial) = nan;
                    end
                end
            end
            cr_habit_remap = data_cr_symb_old(data_cr_symb_old.is_remapped ~= 1,:);
            cr_habit_remap_trials = cr_habit_remap(cr_habit_remap.correct == 1,:);
            Data_online(g).CR.rt_unchange{s,1} = nanmean(cr_habit_remap_trials.rt);
            Data_online(g).CR.N_unchange{s,1} = length(cr_habit_remap_trials.rt);
            
            if ~isempty(data_cr_symb_remap)
                for trial = 1:size(data_cr_symb_remap,1)
                    if data_cr_symb_remap.proposed_choice(trial) == remap_pair(1,1)
                        data_cr_symb_remap.is_remapped(trial) = 1;
                        data_cr_symb_remap.remapped_from(trial) = remap_pair(1,2);
                    elseif data_cr_symb_remap.proposed_choice(trial) == remap_pair(1,2)
                        data_cr_symb_remap.is_remapped(trial) = 1;
                        data_cr_symb_remap.remapped_from(trial) = remap_pair(1,1);
                    elseif data_cr_symb_remap.proposed_choice(trial) == remap_pair(2,1)
                        data_cr_symb_remap.is_remapped(trial) = 1;
                        data_cr_symb_remap.remapped_from(trial) = remap_pair(2,2);
                    elseif data_cr_symb_remap.proposed_choice(trial) == remap_pair(2,2)
                        data_cr_symb_remap.is_remapped(trial) = 1;
                        data_cr_symb_remap.remapped_from(trial) = remap_pair(2,1);
                    else
                        data_cr_symb_remap.is_remapped(trial) = 0;
                        data_cr_symb_remap.remapped_from(trial) = nan;
                    end
                end
             end
            cr_habit_remap = data_cr_symb_remap(data_cr_symb_remap.is_remapped ~= 1,:);
            cr_habit_remap_trials = cr_habit_remap(cr_habit_remap.correct == 1,:);
            Data_online(g).CR.rt_unchange{s,2} = nanmean(cr_habit_remap_trials.rt);
            Data_online(g).CR.N_unchange{s,2} = length(cr_habit_remap_trials.rt);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% plot Phit for the testing session %%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Phit for habit test: forced-RT        
    
        index_blk_type = (all(char(data.block_type) == 'TR',2));
        index_stim_type = (all(char(data.stim_type) == 'Symb',2));

        index_no_change = (data.remap == 0);
        index_new = (data.remap == 1);
        index_remap = (data.grp_swap == 1);
        index_block = (data.block_num < 20);
        data_tr_lett_prob = [];
        
        data_tr_lett_prob_remap = data(index_blk_type == 1 & index_stim_type == 1 & index_new == 1 & index_remap == 1 & index_block == 1,:);

        % deal with bad trials; no response
        data_tr_lett_prob_remap(data_tr_lett_prob_remap.actual_choice == num_trg,:) = [];
        % caluclate true allowed RT
        data_tr_lett_prob_remap.real_prep_time = data_tr_lett_prob_remap.prep_time + data_tr_lett_prob_remap.rt - max_rt;
            
       % creat index for element whether mapped or remapped
        if ~isempty(data_tr_lett_prob_remap)
            for trial = 1:size(data_tr_lett_prob_remap,1)
                if data_tr_lett_prob_remap.proposed_choice(trial) == remap_pair(1,1)
                    data_tr_lett_prob_remap.is_remapped(trial) = 1;
                    data_tr_lett_prob_remap.remapped_from(trial) = remap_pair(1,2);
                elseif data_tr_lett_prob_remap.proposed_choice(trial) == remap_pair(1,2)
                    data_tr_lett_prob_remap.is_remapped(trial) = 1;
                    data_tr_lett_prob_remap.remapped_from(trial) = remap_pair(1,1);
                elseif data_tr_lett_prob_remap.proposed_choice(trial) == remap_pair(2,1)
                    data_tr_lett_prob_remap.is_remapped(trial) = 1;
                    data_tr_lett_prob_remap.remapped_from(trial) = remap_pair(2,2);
                elseif data_tr_lett_prob_remap.proposed_choice(trial) == remap_pair(2,2)
                    data_tr_lett_prob_remap.is_remapped(trial) = 1;
                    data_tr_lett_prob_remap.remapped_from(trial) = remap_pair(2,1);
                else
                    data_tr_lett_prob_remap.is_remapped(trial) = 0;
                    data_tr_lett_prob_remap.remapped_from(trial) = nan;
                end
            end
            % trials with non-remapped stimuli
            unchange_remap = data_tr_lett_prob_remap(data_tr_lett_prob_remap.is_remapped == 0,:); % all unchanged trials
            % trials with remapped stimuli
            remap_remap = data_tr_lett_prob_remap(data_tr_lett_prob_remap.is_remapped == 1,:); % all remapped trials
            
            % for remap_remap trials
            % habit error
            habit_remap = remap_remap; habit_remap.correct(:) = 0;
            habit_remap.correct(habit_remap.actual_choice == habit_remap.remapped_from | habit_remap.actual_choice + 4 == habit_remap.remapped_from) = 1;

            % other errors
            other_remap = remap_remap; other_remap.correct(:) = 1;
            other_remap.correct(habit_remap.correct == 1 | remap_remap.correct == 1) = 0;
        end

        % save remapped data
        % a colume for error type
        remap_remap.errorType = zeros(size(remap_remap,1),1);
        remap_remap.errorType(remap_remap.correct == 1) = 1;
        remap_remap.errorType(habit_remap.correct == 1) = 2;
        remap_remap.errorType(other_remap.correct == 1) = 3;

        responsedata(g).remap{s,1} = remap_remap; % for modeling; not useful for this paper
        responsedata(g).unchange{s,1} = unchange_remap;

        % SAT by smooth window for each type of responses
        [f N] = sliding_window(unchange_remap.real_prep_time*1000, unchange_remap.correct, xplot,t_size);
        p_unchange = f; clear f;

        [f N] = sliding_window(remap_remap.real_prep_time*1000, remap_remap.correct, xplot,t_size);
        p_remap = f; clear f;

        [f N] = sliding_window(other_remap.real_prep_time*1000, other_remap.correct, xplot,t_size);
        p_other = f/6; p_other_total = f; clear f;

        [f N] = sliding_window(habit_remap.real_prep_time*1000, habit_remap.correct, xplot,t_size);
        p_habit = f; clear f;

       temp = [p_unchange p_remap p_other p_habit]; 
       Data_online(g).Phit{s,1} = temp; % save data for plot group mean later
       clear temp;
       
       % fit the SAT for trials with non-remapped stimuli to get t_min
       [para, ycdf] = fit_unchanged(unchange_remap.real_prep_time, unchange_remap.correct);
        Data_online(g).Para{s,1} = para;
        Data_online(g).Ycdf{s,1} = ycdf;
            
        chance = para(3);
        threshold = chance + 0.05;
        ind = find(ycdf >= threshold, 1); %$ t_min

        permutation_interval = [323, 602] - 300; % determined by permutation test; roughly 0 to 300 ms after t_min
            
        % calcualte habit index for the above interval versus interval before t_min   
        interval1 = find(habit_remap.real_prep_time >= 1/1000 & habit_remap.real_prep_time <= ind/1000);
        interval2 = find(habit_remap.real_prep_time >= (ind + 1)/1000 & habit_remap.real_prep_time <= (ind + 300)/1000);

        habit_index_1 = (sum(habit_remap.correct(interval1)))/length(habit_remap.correct(interval1)) - 0.125;
        habit_index_2 = (sum(habit_remap.correct(interval2)))/length(habit_remap.correct(interval2)) - 0.125;

        habit_index = [habit_index_1 habit_index_2];
        Data_online(g).habit_index_chance{s,1} = habit_index; % save data for plot group mean later
        clear habit_index;

        % repeat for each hand
        % dominant hand first
          if unchange_remap.handedness(1) == 'r'
            index_dh = (char(unchange_remap.key) == 'h' |char(unchange_remap.key) == 'u'...
               | char(unchange_remap.key) == 'i'| char(unchange_remap.key) == 'l');
            index_ndh = (char(unchange_remap.key) == 'a' |char(unchange_remap.key) == 'w'...
               | char(unchange_remap.key) == 'e'| char(unchange_remap.key) == 'f');
          
          elseif unchange_remap.handedness(1) == 'l'
            index_dh = (char(unchange_remap.key) == 'a' |char(unchange_remap.key) == 'w'...
               | char(unchange_remap.key) == 'e'| char(unchange_remap.key) == 'f');
            index_ndh = (char(unchange_remap.key) == 'h' |char(unchange_remap.key) == 'u'...
               | char(unchange_remap.key) == 'i'| char(unchange_remap.key) == 'l');
          end
          
          unchange_remap_dh = unchange_remap(index_dh == 1,:);
          unchange_remap_ndh = unchange_remap(index_ndh == 1,:);
          
          % nondominant hand first
          if remap_remap.handedness(1) == 'r'
            index_dh = (char(remap_remap.key) == 'h' |char(remap_remap.key) == 'u'...
               | char(remap_remap.key) == 'i'| char(remap_remap.key) == 'l');
            index_ndh = (char(remap_remap.key) == 'a' |char(remap_remap.key) == 'w'...
               | char(remap_remap.key) == 'e'| char(remap_remap.key) == 'f');
          
          elseif remap_remap.handedness(1) == 'l'
            index_dh = (char(remap_remap.key) == 'a' |char(remap_remap.key) == 'w'...
               | char(remap_remap.key) == 'e'| char(remap_remap.key) == 'f');
            index_ndh = (char(remap_remap.key) == 'h' |char(remap_remap.key) == 'u'...
               | char(remap_remap.key) == 'i'| char(remap_remap.key) == 'l');
          end
          remap_remap_dh = remap_remap(index_dh == 1,:);
          remap_remap_ndh = remap_remap(index_ndh == 1,:);
          
          if other_remap.handedness(1) == 'r'
            index_dh = (char(other_remap.key) == 'h' |char(other_remap.key) == 'u'...
               | char(other_remap.key) == 'i'| char(other_remap.key) == 'l');
            index_ndh = (char(other_remap.key) == 'a' |char(other_remap.key) == 'w'...
               | char(other_remap.key) == 'e'| char(other_remap.key) == 'f');
          
          elseif other_remap.handedness(1) == 'l'
            index_dh = (char(other_remap.key) == 'a' |char(other_remap.key) == 'w'...
               | char(other_remap.key) == 'e'| char(other_remap.key) == 'f');
            index_ndh = (char(other_remap.key) == 'h' |char(other_remap.key) == 'u'...
               | char(other_remap.key) == 'i'| char(other_remap.key) == 'l');
          end
          other_remap_dh = other_remap(index_dh == 1,:);
          other_remap_ndh = other_remap(index_ndh == 1,:);
          
          if habit_remap.handedness(1) == 'r'
            index_dh = (char(habit_remap.key) == 'h' |char(habit_remap.key) == 'u'...
               | char(habit_remap.key) == 'i'| char(habit_remap.key) == 'l');
            index_ndh = (char(habit_remap.key) == 'a' |char(habit_remap.key) == 'w'...
               | char(habit_remap.key) == 'e'| char(habit_remap.key) == 'f');
          
          elseif habit_remap.handedness(1) == 'l'
            index_dh = (char(habit_remap.key) == 'a' |char(habit_remap.key) == 'w'...
               | char(habit_remap.key) == 'e'| char(habit_remap.key) == 'f');
            index_ndh = (char(habit_remap.key) == 'h' |char(habit_remap.key) == 'u'...
               | char(habit_remap.key) == 'i'| char(habit_remap.key) == 'l');
          end
          habit_remap_dh = habit_remap(index_dh == 1,:);
          habit_remap_ndh = habit_remap(index_ndh == 1,:);
          
          
          [f N] = sliding_window(unchange_remap_dh.real_prep_time*1000, unchange_remap_dh.correct, xplot,t_size);
            p_unchange = f; clear f;

            [f N] = sliding_window(remap_remap_dh.real_prep_time*1000, remap_remap_dh.correct, xplot,t_size);
            p_remap = f; clear f;

            [f N] = sliding_window(other_remap_dh.real_prep_time*1000, other_remap_dh.correct, xplot,t_size);
            p_other = f/6; p_other_total = f; clear f;

            [f N] = sliding_window(habit_remap_dh.real_prep_time*1000, habit_remap_dh.correct, xplot,t_size);
            p_habit = f; clear f;

           temp = [p_unchange p_remap p_other p_habit]; 
           Data_online(g).Phit_dh{s,1} = temp; % save data for plot group mean later
           clear temp;
           
           [f N] = sliding_window(unchange_remap_ndh.real_prep_time*1000, unchange_remap_ndh.correct, xplot,t_size);
            p_unchange = f; clear f;

            [f N] = sliding_window(remap_remap_ndh.real_prep_time*1000, remap_remap_ndh.correct, xplot,t_size);
            p_remap = f; clear f;

            [f N] = sliding_window(other_remap_ndh.real_prep_time*1000, other_remap_ndh.correct, xplot,t_size);
            p_other = f/6; p_other_total = f; clear f;

            [f N] = sliding_window(habit_remap_ndh.real_prep_time*1000, habit_remap_ndh.correct, xplot,t_size);
            p_habit = f; clear f;

           temp = [p_unchange p_remap p_other p_habit]; 
           Data_online(g).Phit_ndh{s,1} = temp; % save data for plot group mean later
           clear temp;

%% SAT for set-size effect: forced-RT
            index_blk_type = (all(char(data.block_type) == 'TR',2));
            index_stim_type = (all(char(data.stim_type) == 'Symb',2));

            index_no_change = (data.remap == 0);
            index_new = (data.remap == 1);
            index_remap = (data.grp_swap == 1);
            index_block = (data.block_num < 10);

            ind_set_4 = (data.set_size == 4);
            ind_set_8 = (data.set_size == 8);
        
        ind_stim = (data.proposed_choice == remap_pair(1,1) | data.proposed_choice == remap_pair(1,2) ...
            | data.proposed_choice == remap_pair(2,1) | data.proposed_choice == remap_pair(2,2));
        % find out Tr trials; Rt trials; Cr trials.
        
        data_tr_4 = data(index_blk_type == 1 & index_stim_type == 1 ...
        & index_no_change == 1 & index_block == 1 & ind_set_4 == 1,:);
    
        data_tr_8 = data(index_blk_type == 1 & index_stim_type == 1 ...
        & index_no_change == 1 & index_block == 1 & ind_set_8 & ind_stim == 1,:);
        
        %deal with bad trials
        data_tr_4(data_tr_4.actual_choice == num_trg,:) = [];
        data_tr_8(data_tr_8.actual_choice == num_trg,:) = [];
        
        data_tr_4.real_prep_time = data_tr_4.prep_time + data_tr_4.rt - max_rt;
        data_tr_8.real_prep_time = data_tr_8.prep_time + data_tr_8.rt - max_rt;
        
        % SAT smooth window
        [f N] = sliding_window(data_tr_4.real_prep_time*1000, data_tr_4.correct, xplot,t_size);
        p_tr_4 = f; clear f;
        
        [f N] = sliding_window(data_tr_8.real_prep_time*1000, data_tr_8.correct, xplot,t_size);
        p_tr_8 = f; clear f;

        temp = [p_tr_4 p_tr_8];
        Data_online(g).SAT{s,1} = temp; % save data for plot group mean later
        clear temp;
        
        % fit to get mean speed \mu
        [para, ycdf] = fit_unchanged(data_tr_4.real_prep_time, data_tr_4.correct);
        Data_online(g).SAT_Para{s,1} = para;
        Data_online(g).SAT_Ycdf{s,1} = ycdf;
        
        [para, ycdf] = fit_unchanged(data_tr_8.real_prep_time, data_tr_8.correct);
        Data_online(g).SAT_Para{s,2} = para;
        Data_online(g).SAT_Ycdf{s,2} = ycdf;
        
        % deal with repetition effect
        % remove repetition trials 
        % n and n+1 need the same response
        data_tr_4_sort = sortrows(data_tr_4,{'block_num','trial_Count'});
        data_tr_8_sort = sortrows(data_tr_8,{'block_num','trial_Count'});
        
        index_4_nrp = find(diff(data_tr_4_sort.proposed_choice) ~= 0) + 1;
        index_8_nrp = find(diff(data_tr_8_sort.proposed_choice) ~= 0) + 1;
        
        n_4_nrp = (size(data_tr_4_sort,1) - length(index_4_nrp))/size(data_tr_4_sort,1);
        n_8_nrp = (size(data_tr_8_sort,1) - length(index_8_nrp))/size(data_tr_8_sort,1);
        
        data_tr_4_nrp = data_tr_4_sort(index_4_nrp,:);
        data_tr_8_nrp = data_tr_8_sort(index_8_nrp,:);
        
        [f N] = sliding_window(data_tr_4_nrp.real_prep_time*1000, data_tr_4_nrp.correct, xplot,t_size);
        p_tr_4 = f; clear f;
        
        [f N] = sliding_window(data_tr_8_nrp.real_prep_time*1000, data_tr_8_nrp.correct, xplot,t_size);
        p_tr_8 = f; clear f;

        temp = [p_tr_4 p_tr_8];
        Data_online(g).SAT_nrp{s,1} = temp; % save data for plot group mean later
        clear temp;
        
        [para, ycdf] = fit_unchanged(data_tr_4_nrp.real_prep_time, data_tr_4_nrp.correct);
        Data_online(g).SAT_Para_nrp{s,1} = para;
        Data_online(g).SAT_Ycdf_nrp{s,1} = ycdf;
        
        [para, ycdf] = fit_unchanged(data_tr_8_nrp.real_prep_time, data_tr_8_nrp.correct);
        Data_online(g).SAT_Para_nrp{s,2} = para;
        Data_online(g).SAT_Ycdf_nrp{s,2} = ycdf;
        
        Data_online(g).SAT_nrp_n{s,1} = n_4_nrp;
        Data_online(g).SAT_nrp_n{s,2} = n_8_nrp;
        
        % repeat above process for handedness
        if data_tr_4.handedness(1) == 'r'
            index_dh = (char(data_tr_4.key) == 'h' |char(data_tr_4.key) == 'u'...
               | char(data_tr_4.key) == 'i'| char(data_tr_4.key) == 'l');
            index_ndh = (char(data_tr_4.key) == 'a' |char(data_tr_4.key) == 'w'...
               | char(data_tr_4.key) == 'e'| char(data_tr_4.key) == 'f');
          
          elseif data_tr_4.handedness(1) == 'l'
            index_dh = (char(data_tr_4.key) == 'a' |char(data_tr_4.key) == 'w'...
               | char(data_tr_4.key) == 'e'| char(data_tr_4.key) == 'f');
            index_ndh = (char(data_tr_4.key) == 'h' |char(data_tr_4.key) == 'u'...
               | char(data_tr_4.key) == 'i'| char(data_tr_4.key) == 'l');
        end
          
        data_tr_4_dh = data_tr_4(index_dh == 1,:);
        data_tr_4_ndh = data_tr_4(index_ndh == 1,:);
        
        clear index_dh index_ndh;
        if data_tr_8.handedness(1) == 'r'
            index_dh = (char(data_tr_8.key) == 'h' |char(data_tr_8.key) == 'u'...
               | char(data_tr_8.key) == 'i'| char(data_tr_8.key) == 'l');
            index_ndh = (char(data_tr_8.key) == 'a' |char(data_tr_8.key) == 'w'...
               | char(data_tr_8.key) == 'e'| char(data_tr_8.key) == 'f');
          
          elseif data_tr_8.handedness(1) == 'l'
            index_dh = (char(data_tr_8.key) == 'a' |char(data_tr_8.key) == 'w'...
               | char(data_tr_8.key) == 'e'| char(data_tr_8.key) == 'f');
            index_ndh = (char(data_tr_8.key) == 'h' |char(data_tr_8.key) == 'u'...
               | char(data_tr_8.key) == 'i'| char(data_tr_8.key) == 'l');
        end
        data_tr_8_dh = data_tr_8(index_dh == 1,:);
        data_tr_8_ndh = data_tr_8(index_ndh == 1,:);
        
        [f N] = sliding_window(data_tr_4_dh.real_prep_time*1000, data_tr_4_dh.correct, xplot,t_size);
        p_tr_4 = f; clear f;
        
        [f N] = sliding_window(data_tr_8_dh.real_prep_time*1000, data_tr_8_dh.correct, xplot,t_size);
        p_tr_8 = f; clear f;

        temp = [p_tr_4 p_tr_8];
        Data_online(g).SAT_dh{s,1} = temp; % save data for plot group mean later
        clear temp;
        
        [para, ycdf] = fit_unchanged(data_tr_4_dh.real_prep_time, data_tr_4_dh.correct);
        Data_online(g).SAT_Para_dh{s,1} = para;
        Data_online(g).SAT_Ycdf_dh{s,1} = ycdf;
        
        [para, ycdf] = fit_unchanged(data_tr_8_dh.real_prep_time, data_tr_8_dh.correct);
        Data_online(g).SAT_Para_dh{s,2} = para;
        Data_online(g).SAT_Ycdf_dh{s,2} = ycdf;
        
        [f N] = sliding_window(data_tr_4_ndh.real_prep_time*1000, data_tr_4_ndh.correct, xplot,t_size);
        p_tr_4 = f; clear f;
        
        [f N] = sliding_window(data_tr_8_ndh.real_prep_time*1000, data_tr_8_ndh.correct, xplot,t_size);
        p_tr_8 = f; clear f;

        temp = [p_tr_4 p_tr_8];
        Data_online(g).SAT_ndh{s,1} = temp; % save data for plot group mean later
        clear temp;
        
        [para, ycdf] = fit_unchanged(data_tr_4_ndh.real_prep_time, data_tr_4_ndh.correct);
        Data_online(g).SAT_Para_ndh{s,1} = para;
        Data_online(g).SAT_Ycdf_ndh{s,1} = ycdf;
        
        [para, ycdf] = fit_unchanged(data_tr_8_ndh.real_prep_time, data_tr_8_ndh.correct);
        Data_online(g).SAT_Para_ndh{s,2} = para;
        Data_online(g).SAT_Ycdf_ndh{s,2} = ycdf;
        
%% set-size effect; free-RT
        ind_blk_type = (all(char(data.block_type) == 'RT',2));
        ind_stim_type = (all(char(data.stim_type) == 'Symb',2));
        ind_pre_rt = (data.pre_rt == 1);
        ind_set_4 = (data.set_size == 4);
        ind_set_8 = (data.set_size == 8);
        
        ind_stim = (data.proposed_choice == remap_pair(1,1) | data.proposed_choice == remap_pair(1,2) ...
            | data.proposed_choice == remap_pair(2,1) | data.proposed_choice == remap_pair(2,2));
        data_rt_symb_4 = data(ind_blk_type == 1 & ind_stim_type == 1 & ind_pre_rt == 0 & ind_set_4 == 1,:);
        data_rt_symb_8 = data(ind_blk_type == 1 & ind_stim_type == 1 & ind_pre_rt == 0 & ind_set_8 == 1 & ind_stim == 1,:);
        
        data_rt_symb_4(data_rt_symb_4.rt > 3,:) = [];
        data_rt_symb_8(data_rt_symb_8.rt > 3,:) = [];

        mean_prep_time_4 = nanmean(data_rt_symb_4.rt(data_rt_symb_4.correct == 1));
        mean_correct_4 = nanmean(data_rt_symb_4.correct);

        mean_prep_time_8 = nanmean(data_rt_symb_8.rt(data_rt_symb_8.correct == 1));
        mean_correct_8 = nanmean(data_rt_symb_8.correct);
        
        Data_online(g).RT{s,1} = mean_prep_time_4;
        Data_online(g).RT{s,2} = mean_prep_time_8;
        
        Data_online(g).AC{s,1} = mean_correct_4;
        Data_online(g).AC{s,2} = mean_correct_8;
        
        % no repetition trials
        data_rt_symb_4_sort = sortrows(data_rt_symb_4,{'block_num','trial_Count'});
        data_rt_symb_8_sort = sortrows(data_rt_symb_8,{'block_num','trial_Count'});
        
        index_4_nrp = find(diff(data_rt_symb_4_sort.proposed_choice) ~= 0) + 1;
        index_8_nrp = find(diff(data_rt_symb_8_sort.proposed_choice) ~= 0) + 1;
        
        n_4_nrp = (size(data_rt_symb_4_sort,1) - length(index_4_nrp))/size(data_rt_symb_4_sort,1);
        n_8_nrp = (size(data_rt_symb_8_sort,1) - length(index_8_nrp))/size(data_rt_symb_8_sort,1);
        
        data_rt_symb_4_nrp = data_rt_symb_4_sort(index_4_nrp,:);
        data_rt_symb_8_nrp = data_rt_symb_8_sort(index_8_nrp,:);
        
        mean_prep_time_4 = nanmean(data_rt_symb_4_nrp.rt(data_rt_symb_4_nrp.correct == 1));
        mean_correct_4 = nanmean(data_rt_symb_4_nrp.correct);

        mean_prep_time_8 = nanmean(data_rt_symb_8_nrp.rt(data_rt_symb_8_nrp.correct == 1));
        mean_correct_8 = nanmean(data_rt_symb_8_nrp.correct);
        
        Data_online(g).RT_nrp{s,1} = mean_prep_time_4;
        Data_online(g).RT_nrp{s,2} = mean_prep_time_8;
        
        Data_online(g).AC_nrp{s,1} = mean_correct_4;
        Data_online(g).AC_nrp{s,2} = mean_correct_8;
        
        Data_online(g).RT_nrp_n{s,1} = n_4_nrp;
        Data_online(g).RT_nrp_n{s,2} = n_8_nrp;
        
        % repeat for handedness
        if data_rt_symb_4.handedness(1) == 'r'
            index_dh = (char(data_rt_symb_4.key) == 'h' |char(data_rt_symb_4.key) == 'u'...
               | char(data_rt_symb_4.key) == 'i'| char(data_rt_symb_4.key) == 'l');
            index_ndh = (char(data_rt_symb_4.key) == 'a' |char(data_rt_symb_4.key) == 'w'...
               | char(data_rt_symb_4.key) == 'e'| char(data_rt_symb_4.key) == 'f');
          
          elseif data_rt_symb_4.handedness(1) == 'l'
            index_dh = (char(data_rt_symb_4.key) == 'a' |char(data_rt_symb_4.key) == 'w'...
               | char(data_rt_symb_4.key) == 'e'| char(data_rt_symb_4.key) == 'f');
            index_ndh = (char(data_rt_symb_4.key) == 'h' |char(data_rt_symb_4.key) == 'u'...
               | char(data_rt_symb_4.key) == 'i'| char(data_rt_symb_4.key) == 'l');
        end
        data_rt_symb_4_dh = data_rt_symb_4(index_dh == 1,:);
        data_rt_symb_4_ndh = data_rt_symb_4(index_ndh == 1,:);
        
        if data_rt_symb_8.handedness(1) == 'r'
            index_dh = (char(data_rt_symb_8.key) == 'h' |char(data_rt_symb_8.key) == 'u'...
               | char(data_rt_symb_8.key) == 'i'| char(data_rt_symb_8.key) == 'l');
            index_ndh = (char(data_rt_symb_8.key) == 'a' |char(data_rt_symb_8.key) == 'w'...
               | char(data_rt_symb_8.key) == 'e'| char(data_rt_symb_8.key) == 'f');
          
          elseif data_rt_symb_8.handedness(1) == 'l'
            index_dh = (char(data_rt_symb_8.key) == 'a' |char(data_rt_symb_8.key) == 'w'...
               | char(data_rt_symb_8.key) == 'e'| char(data_rt_symb_8.key) == 'f');
            index_ndh = (char(data_rt_symb_8.key) == 'h' |char(data_rt_symb_8.key) == 'u'...
               | char(data_rt_symb_8.key) == 'i'| char(data_rt_symb_8.key) == 'l');
        end
        data_rt_symb_8_dh = data_rt_symb_8(index_dh == 1,:);
        data_rt_symb_8_ndh = data_rt_symb_8(index_ndh == 1,:);
        
        mean_prep_time_4 = nanmean(data_rt_symb_4_dh.rt(data_rt_symb_4_dh.correct == 1));
        mean_correct_4 = nanmean(data_rt_symb_4_dh.correct);

        mean_prep_time_8 = nanmean(data_rt_symb_8_dh.rt(data_rt_symb_8_dh.correct == 1));
        mean_correct_8 = nanmean(data_rt_symb_8_dh.correct);
        
        Data_online(g).RT_dh{s,1} = mean_prep_time_4;
        Data_online(g).RT_dh{s,2} = mean_prep_time_8;
        
        Data_online(g).AC_dh{s,1} = mean_correct_4;
        Data_online(g).AC_dh{s,2} = mean_correct_8;
        
        mean_prep_time_4 = nanmean(data_rt_symb_4_ndh.rt(data_rt_symb_4_ndh.correct == 1));
        mean_correct_4 = nanmean(data_rt_symb_4_ndh.correct);

        mean_prep_time_8 = nanmean(data_rt_symb_8_ndh.rt(data_rt_symb_8_ndh.correct == 1));
        mean_correct_8 = nanmean(data_rt_symb_8_ndh.correct);
        
        Data_online(g).RT_ndh{s,1} = mean_prep_time_4;
        Data_online(g).RT_ndh{s,2} = mean_prep_time_8;
        
        Data_online(g).AC_ndh{s,1} = mean_correct_4;
        Data_online(g).AC_ndh{s,2} = mean_correct_8;
    end
end
clear DATA;

DATA_NA.Online = Data_online;
datafname = ['Data_for_Analysis'];
save(datafname, 'DATA_NA');

RESPONSE.Online = responsedata;
datafname = ['Data_for_Model'];
save(datafname,'RESPONSE');