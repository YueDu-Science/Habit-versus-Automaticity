function [C_global] = permutation_test_2tailed(x,y, paired, grp)

% reference1: Safaie et al., 2019 Turning the body into a clock: Accurate timing is
%facilitated by simple stereotyped interactions
% with the environment - SI

%reference 2: Fujisawa et al., 2008 Behavior-dependent short-term assembly dynamics in
% the medial prefrontal cortex

% grp: num of groups: adjust alpha level so that family-wise type I error rate is 0.05;

% x is habitual error
% y is other errors
num = size(x,2); 
if paired == 1
    D(1,:) = nanmean(x' - y');
else
    D(1,:) = nanmean(x') - nanmean(y');
end

alpha = 0.05/grp;

% first step
% permute data: assign group lables to each data set
% 0: habitual errors
% 1: other errors

% since our data are paired within each individual
% we randomly assign group lable to other errors and habitual errors
n = 10001;
for i = 2:n
    grp = randi(2, num, 1) - 1; % generate random 0 and 1
    ind0 = find(grp == 0); ind1 = find(grp == 1);
    x_perm = [x(:,ind0) y(:,ind1)];
    y_perm = [y(:,ind0) x(:,ind1)];
    if paired == 1
        D(i,:) = nanmean(x_perm' - y_perm');
    elseif paired == 0
        D(i,:) = nanmean(x_perm') - nanmean(y_perm');
    end
end

% calculate point wise p-value
T = size(x,1);
% for t = 1:T
%     pv_point_plus(t) = sum(D(2:end,t) >= D(1,t))/numel(D(:,t));
%     pv_point_minus(t) = sum(D(2:end,t) <= D(1,t))/numel(D(:,t));
%     pv_point(t) = min([1,2*pv_point_plus, 2*pv_point_minus]);
% end

% point-wise band
c_range = [min(D); max(D)];
% for t = 1:T
%     step_range = linspace(0, c_range(2,t), 200);
%     for m = 1:length(step_range)
%         if (sum(D(:,t) >= step_range(m))/numel(D(:,t))) <= alpha/2 %|| (sum(D(:,t) <= -step_range(m))/numel(D(:,t))) <= alpha
%             continue;
%         else
%             C_point_upper(t) = step_range(m); 
%         end
%     end
% end

% for t = 1:T
%     step_range = linspace(0, c_range(1,t), 200);
%     for m = 1:length(step_range)
%         if (sum(D(:,t) <= step_range(m))/numel(D(:,t))) <= alpha/2 %|| (sum(D(:,t) <= -step_range(m))/numel(D(:,t))) <= alpha
%             continue;
%         else
%             C_point_lower(t) = step_range(m); 
%         end
%     end
% end
% C_point = [C_point_upper; C_point_lower];

% to generate multiple point wise band with different alpha = a
 % 0.05 is the alpha level without correction to multiple comparison
% since we are using two-tailed test (whether habitual > other errors)
% we need to find c such than the frequency of that data is larger than c
% is smaller than a
% for a two-tailed test, find cc that the frequency of that data is
% larger than a/2 and smaller than -c is smaller than a/2

% in other words, construct a  line that at each response time
% point, the value is larger than 95% of the permuted data

a = linspace(0.01, 0.00001, 3);  %if stopped after iter 0; then increase the first value.
success = 0;
iter = 0
while success == 0
    C_upper = [];
    for j = 1:length(a)
        for t = 1:T
            step_range = linspace(0, c_range(2,t), 200);
            for m = 1:length(step_range)
                if (sum(D(:,t) >= step_range(m))/numel(D(:,t))) <= a(j)/2 %|| (sum(D(:,t) <= -step_range(m))/numel(D(:,t))) > a(j)
                    continue;
                else
                    C_upper(j,t) = step_range(m); 
                end
            end
        end
    end
    
    C_lower = [];
    for j = 1:length(a)
        for t = 1:T
            step_range = linspace(0, c_range(1,t), 200);
            for m = 1:length(step_range)
                if (sum(D(:,t) <= step_range(m))/numel(D(:,t))) <= a(j)/2 %|| (sum(D(:,t) <= -step_range(m))/numel(D(:,t))) > a(j)
                    continue;
                else
                    C_lower(j,t) = step_range(m); 
                end
            end
        end
    end
    %C(1,:) is the point wise 0.05 level

    % to determine which C(j,:) has family wise error as alpha = 0.05
    count_upper = [];
    count_lower = [];
    for j = 1:size(C_upper,1)
        for n = 1:size(D,1)
            if sum(D(n,:) > C_upper(j,:)) > 0   % if across all time point, there data points beyond the c_upper or c_lower
                count_upper(n,j) = 1;
            else
                count_upper(n,j) = 0;
            end
            if sum(D(n,:) < C_lower(j,:)) > 0
                count_lower(n,j) = 1;
            else
                count_lower(n,j) = 0;
            end
        end
    end
    prob_upper = nanmean(count_upper);
    prob_lower = nanmean(count_lower);
    ind_upper = find(prob_upper < alpha/2,1);
    ind_lower = find(prob_lower < alpha/2,1);
    
    if ~isempty(ind_upper) && isempty(ind_lower)
       ind = ind_upper;
    elseif isempty(ind_upper) && ~isempty(ind_lower)
        ind = ind_lower;
    elseif ~isempty(ind_upper) && ~isempty(ind_lower)
        ind = min(ind_upper, ind_lower);
    end
    
    if ind == 1
        break;
    else
        if a(ind-1) - a(ind) < 1e-5
            success = 1;
        else
            success = 0;
            a = linspace(a(ind-1), a(ind), 3);
        end
        iter = iter + 1
    end
    
end

C_global_upper = C_upper(ind,:);
C_global_lower = C_lower(ind,:);

C_global = [C_global_upper; C_global_lower];

% plot(D(1,:),'k-')
% plot(C_point,'b--');
% plot(C(ind,:),'r--')

%count(ind,:) is the family wise alpha level