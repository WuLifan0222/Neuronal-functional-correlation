function [C_select_by_trials_ind,C_select_by_trials,C_trials,neuron_trial_similar,all_neuron_response_for_each_sti_condition] = select_by_trials(C,param,stimuli)
sti_num = param.sti_num;
trial_threshold = param.trial_select_thresh;
start_edge = stimuli.start_edge;
end_edge = stimuli.end_edge;

neuron_trial_similar = []; 
C_trials = {};
all_neuron_response_bool = [];
all_neuron_response_for_each_sti_condition = [];
for i = 1:size(C,1)
    neuron_single = C(i,:);
    [trials_cell, before, during, after, sti_labels,sti_num, trials_num,response_bool,response_for_each_sti_condition] = ...
        get_multi_neurons_trials_matrix_by_sti_with_sel(neuron_single, stimuli, param);
    all_neuron_response_bool = [all_neuron_response_bool,response_bool];
    all_neuron_response_for_each_sti_condition = [all_neuron_response_for_each_sti_condition;response_for_each_sti_condition];
    C_trials{i} = trials_cell;
end
select_first_id = find(all_neuron_response_bool==1);
% % % % C_select_first = C(select_first_id,:);

neuron_trial_similar = []; 
for ii = 1:length(select_first_id)
    trials_cell = C_trials{ii};
    for sti_i = 1:sti_num
        similar = [];
        similar_matrix = [];
        C_trials_matrix = trials_cell{sti_i};
        Y = corrcoef(C_trials_matrix);
        similar_for_single_sti = sum(sum(Y));
        similar_matrix(sti_i) = similar_for_single_sti;
        similar = mean(similar_matrix);
    end
    neuron_trial_similar = [neuron_trial_similar,similar];
end

neuron_trial_similar_nonan = neuron_trial_similar(~isnan(neuron_trial_similar));
select_second_id = select_first_id(find(~isnan(neuron_trial_similar)));
C_select_second = C(select_second_id,:);

[neuron_sort,neuron_sort_ind] = sort(neuron_trial_similar_nonan,'descend');
C_select_by_trials_ind = select_second_id(neuron_sort_ind);
C_select_by_trials = C(C_select_by_trials_ind,:);
