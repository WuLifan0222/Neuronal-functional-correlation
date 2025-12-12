function [edge_stdResult,sort_edge_stdResult_id,sort_edge_stdResult_value]= get_all_stdResult_value(select_neurons_C_sort_by_region_corr_matrix,stimuli, param)
trial_frames = (param.sti_during_time + param.show_sti_after) * param.fs;
[select_neurons_C_sort_by_region_corr_matrix_cell,stimuli_sort,select_neurons_C_sort_by_region_corr_matrix_visualsort,~, ~, ~, ~,~, ~] = ...
    sort_multi_neurons_trials_matrix_by_sti(select_neurons_C_sort_by_region_corr_matrix, stimuli, param); 
edge_stdResult = {};
sort_edge_stdResult_id = {};
sort_edge_stdResult_value = {};

for index0 = 1:length(select_neurons_C_sort_by_region_corr_matrix_cell)
    nNeuron = size(select_neurons_C_sort_by_region_corr_matrix_cell{index0},2);
    ori_select_neurons_C_R = reshape(select_neurons_C_sort_by_region_corr_matrix_cell{index0}', nNeuron, trial_frames, []);
    % ori_select_neurons_C_R = reshape(select_neurons_C_sort_by_region_corr_matrix_cell{index0}, nNeuron, trial_frames, []);
    stdResult_buff = zeros(nNeuron,1); 
    for neuroni = 1:nNeuron
        % stdResult(i) = calResNormStd(ori_select_neurons_C_R(i,:, 1:8));
        stdResult_buff(neuroni) = calResNormStd(ori_select_neurons_C_R(neuroni,:, :));    
    end
    [sortedValues_buff, sortedIndices_buff] = sort(stdResult_buff);
    edge_stdResult{index0} = stdResult_buff;
    sort_edge_stdResult_id{index0} = sortedIndices_buff;
    sort_edge_stdResult_value{index0} = sortedValues_buff;
end
