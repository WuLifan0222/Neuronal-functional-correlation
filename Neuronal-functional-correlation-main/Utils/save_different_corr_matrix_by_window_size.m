function corr_param = save_different_corr_matrix_by_window_size(C,select_neurons_C,win_size,param,stimuli,result_for_response_to_stis_path)
corr_param.corr_win = win_size;
corr_param.corr_win_before_num = param.show_sti_before*param.fs/corr_param.corr_win;
corr_param.corr_win_during_num = param.sti_during_time*param.fs/corr_param.corr_win;
corr_param.corr_win_after_num = param.show_sti_after*param.fs/corr_param.corr_win;
corr_param.corr_win_num = size(select_neurons_C,2)-corr_param.corr_win;
[select_neurons_C_sort_by_region_corr_matrix, select_neurons_C_sort_by_region_edge_cell] = get_corr_matrix_from_neuron_matrix(select_neurons_C,corr_param);

save([result_for_response_to_stis_path sprintf('neuron_region_lr_with_select_neurons_win%d.mat',win_size)],...
    'select_neurons_C_sort_by_region_corr_matrix',...
    'select_neurons_C_sort_by_region_edge_cell', ...
    'param','stimuli','-v7.3');
clear select_neurons_C_sort_by_region_corr_matrix;
clear select_neurons_C_sort_by_region_edge_cell;
