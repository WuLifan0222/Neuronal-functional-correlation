function [trials_cell, before, during, after, sti_labels,sti_num, trials_num,response_bool,response_for_each_sti_condition] = get_multi_neurons_trials_matrix_by_sti_with_sel(single_trace, stimuli, param)
single_trace = single_trace';
before = round(param.show_sti_before * param.fs);
during = round(param.sti_during_time * param.fs);
after = round(param.show_sti_after * param.fs);
start_edge = stimuli.start_edge;
end_edge = stimuli.end_edge;
sti_labels = unique(stimuli.stimuili_label_ind);
sti_num = length(sti_labels);  
trials_cell = {};  
trials_num = [];
response_for_each_sti_condition = [];
for j = 1:sti_num
    sti_i = sti_labels(j);
    trial_matrix = [];
    trial_before_matrix = [];
    trial_during_matrix = [];
    for i = 1:size(start_edge,2)
        if stimuli.stimuili_label_ind(i) == sti_i
            % trial = single_trace(round(start_edge(i) - before)+1: round(start_edge(i) + during + after), :);
            trial = single_trace(round(start_edge(i) - before): round(start_edge(i) + during + after-1), :);
            trial_matrix = [trial_matrix,trial];
            trial_before = single_trace(round(start_edge(i)- before): round(start_edge(i)-1), :);
            trial_before_matrix = [trial_before_matrix,trial_before];
            trial_during = single_trace(round(start_edge(i)): round(start_edge(i) + during -1), :);
            trial_during_matrix = [trial_during_matrix,trial_during];
        end
    end
    trials_num(j) = size(trial_matrix,2);
    trials_cell{j} = trial_matrix ;

    response_for_each_trial_condition = find(mean(trial_during_matrix)>param.trial_select_thresh*mean(trial_before_matrix));

    if size(response_for_each_trial_condition,2)>=trials_num*param.response_ratio
        response_each_trial_bool = 1;
    else
        response_each_trial_bool = 0;
    end
    response_for_each_sti_condition = [response_for_each_sti_condition,response_each_trial_bool];
end
if sum(response_for_each_sti_condition)>=1
    response_bool = 1;
else
    response_bool = 0;
end

