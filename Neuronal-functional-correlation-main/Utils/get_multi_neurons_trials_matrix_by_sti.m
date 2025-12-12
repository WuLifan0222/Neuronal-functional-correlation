function [trials_cell, before, during, after, sti_labels,sti_num, trials_num] = get_multi_neurons_trials_matrix_by_sti(trace, stimuli, param)
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
for j = 1:sti_num
    sti_i = sti_labels(j);
    trial_matrix = [];
    for i = 1:size(start_edge,2)
        if stimuli.stimuili_label_ind(i) == sti_i
            trial = trace(round(start_edge(i) - before): round(start_edge(i) + during + after-1), :);
            trial_matrix = [trial_matrix,trial];
        end
    end
    trials_num(j) = size(trial_matrix,2);
    trials_cell{j} = trial_matrix ;
end
