function [trials_cell, stimuli_sort,sort_trace,before, during, after, sti_labels,sti_num, trials_num] = sort_multi_neurons_trials_matrix_by_sti(trace, stimuli, param)
trace = trace';
before = round(param.show_sti_before * param.fs);
during = round(param.sti_during_time * param.fs);
after = round(param.show_sti_after * param.fs);
start_edge = stimuli.start_edge;
end_edge = stimuli.end_edge;
sti_labels = unique(stimuli.stimuili_label_ind);
sti_num = length(sti_labels);  
trials_cell = {};  
trials_num = [];
sort_trace = [];
sort_sti = [];
sort_sti_label = [];
sort_sti_label_for_all_trial = [];
label_ind = [];
for j = 1:sti_num
    sti_i = sti_labels(j);
    trial_matrix = [];
    for i = 1:size(start_edge,2)
        if stimuli.stimuili_label_ind(i) == sti_i
            trial = trace(round(start_edge(i) - before): round(start_edge(i) + during + after - before-1), :);
            trial_matrix = [trial_matrix;trial];
            curr_sti_array = [zeros(1,before),ones(1,during),zeros(1,after-before)];
            sort_sti = [sort_sti,curr_sti_array];
            label_ind = [label_ind,sti_i];
        end
    end
    trials_num(j) = size(trial_matrix,2);
    trials_cell{j} = trial_matrix ;
    sort_trace = [sort_trace;trial_matrix];

    sort_sti_label_for_all_trial = [sort_sti_label_for_all_trial,ones(1,size(trial_matrix,1))*j];
    sort_sti_label = sort_sti.*sort_sti_label_for_all_trial;
end

edge = diff(sort_sti_label); 
edge_start = find(edge > 0);
edge_end = edge_start + during;
stimuli_sort.start_edge = edge_start;
stimuli_sort.end_edge = edge_end;
stimuli_sort.stimuili_label = sort_sti_label;
stimuli_sort.stimuili_label_ind = label_ind;
stimuli_sort.stimuli_array = sort_sti;
stimuli_sort.stimuli_array_with_label = sort_sti_label;
stimuli_sort.sort_sti_label_for_all_trial = sort_sti_label_for_all_trial ;
sort_trace = sort_trace';
