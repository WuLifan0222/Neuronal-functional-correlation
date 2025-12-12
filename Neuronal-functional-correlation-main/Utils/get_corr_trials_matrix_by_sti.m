function [corr_trials_cell,before,during,after,sti_num] = get_corr_trials_matrix_by_sti(trace,stimuli,param)
before = param.show_sti_before * param.fs;
during = param.sti_during_time * param.fs;
after = param.show_sti_after * param.fs;
start_edge = stimuli.start_edge;
end_edge = stimuli.end_edge;
sti_labels = unique(stimuli.stimuili_label_ind,'sorted');
sti_num = length(sti_labels);
corr_trials_cell = {};
for j = 1:sti_num
    sti_i = sti_labels(j);
    trials_matrix_for_each_sti = [];
    for i = 1:length(start_edge)          
        if (stimuli.stimuili_label_ind(i)==sti_i)
            trial = [trace(start_edge(i)-before:start_edge(i)+during+after-1)];
            trials_matrix_for_each_sti = [trials_matrix_for_each_sti;trial];   
        end  
    end
    corr_trials_cell{1,sti_i} = trials_matrix_for_each_sti;
end
