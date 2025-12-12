clc; clear; close all; set(0,'defaultfigurecolor','w');
%% analysis nature image sti response  --WLF 20230811
addpath(genpath('D:\PPT_gather\Neuronal functional connectivity\Utils'));

%% >>>>>>>>> load color
color_scheme_npg = get_color();

%% >>>>>>>>> Load experiment (load data and stimuli)
folder_path = ['D:\PPT_gather\Neuronal functional connectivity\dataset\'];
input_C_data = [folder_path 'demo_neurons_data.mat'];
load(input_C_data);
C = C_mat;

param.res_nums_thre = 8; 
result_path_threshold_delay = [folder_path sprintf('result_for_thre%d\\',param.res_nums_thre)];
selected_data = [result_path_threshold_delay 'select_neurons.mat'];
load(selected_data);

% >>>>>>>>> select neuron for all brain region
C_select_for_multi_sti_cell_ind = {};
C_select_for_multi_sti_cell = {};
response_to_stis_neurons = sum(select_bool_for_all_neurons,2);

% >>>>>>>>>>>> statistic neurons response to different stis
for i = 1:param.sti_num
    response_neuron_id = find(select_bool_for_all_neurons(i,:)==1);
    response_neuron_C = C(response_neuron_id,:);
    C_select_for_multi_sti_cell_ind{i} = response_neuron_id;
    C_select_for_multi_sti_cell{i} = response_neuron_C;
end

% >>>>> define select neurons to different regions (act)
select_neurons_id_in_valid_C = C_mat_id_in_validC(select_neurons_id);
region_num = length(brain_region);
select_neurons_sub_left_region_id_invalidC = {};select_neurons_sub_left_region_C = {};
select_neurons_sub_right_region_id_invalidC = {};select_neurons_sub_right_region_C = {};
for i = 1:region_num
    curr_left_brain_id = C_array_left_id{i};
    curr_left_brain_select_id = intersect(curr_left_brain_id, select_neurons_id_in_valid_C);
    select_neurons_sub_left_region_id_invalidC{i} = curr_left_brain_select_id;
    select_neurons_sub_left_region_C{i} = valid_C(curr_left_brain_select_id,:);
    curr_right_brain_id = C_array_right_id{i};
    curr_right_brain_select_id = intersect(curr_right_brain_id, select_neurons_id_in_valid_C);
    select_neurons_sub_right_region_id_invalidC{i} = curr_right_brain_select_id;
    select_neurons_sub_right_region_C{i} = valid_C(curr_right_brain_select_id,:);
end

% >>>>> define select neurons to different regions (deact)
select_deact_neurons_id_in_valid_C = C_mat_id_in_validC(select_deact_neurons_id);
select_deact_neurons_sub_left_region_id_invalidC = {};select_deact_neurons_sub_left_region_C = {};
select_deact_neurons_sub_right_region_id_invalidC = {};select_deact_neurons_sub_right_region_C = {};
for i = 1:region_num
    curr_left_brain_id = C_array_left_id{i};
    curr_left_brain_select_deact_id = intersect(curr_left_brain_id, select_deact_neurons_id_in_valid_C);
    select_deact_neurons_sub_left_region_id_invalidC{i} = curr_left_brain_select_deact_id;
    select_deact_neurons_sub_left_region_C{i} = valid_C(curr_left_brain_select_deact_id,:);
    curr_right_brain_id = C_array_right_id{i};
    curr_right_brain_select_deact_id = intersect(curr_right_brain_id, select_deact_neurons_id_in_valid_C);
    select_deact_neurons_sub_right_region_id_invalidC{i} = curr_right_brain_select_deact_id;
    select_deact_neurons_sub_right_region_C{i} = valid_C(curr_right_brain_select_deact_id,:);
end

save([result_path_threshold_delay 'neuron_region_lr_with_select_neurons.mat'],...
    'valid_C','valid_neuron_x','valid_neuron_y','C_array_id_in_validC','C_mat','C_trials','brain_region_ori','param','stimuli',...
    'unique_brain_region','brain_region','xlabel_cell','neuron_region_id','neuron_num',...
    'neuron_position_x_mat','neuron_position_y_mat','select_neurons_id_in_valid_C','select_deact_neurons_id_in_valid_C',...
    'C_array','xlabel_annoy','neuron_position_x', 'neuron_position_y', ...
    'C_array_left','neuron_position_x_left', 'neuron_position_y_left', 'C_array_left_id',...
    'C_array_right','neuron_position_x_right', 'neuron_position_y_right','C_array_right_id',...
    'select_neurons_C','select_neurons_id','select_deact_neurons_C','select_deact_neurons_id',...
    'C_select_for_multi_sti_cell_ind','C_select_for_multi_sti_cell','select_bool_for_all_neurons','param',...
    'C_select_deact_for_multi_sti_cell_ind','C_select_deact_for_multi_sti_cell','select_deact_bool_for_all_neurons','param',...
    'select_neurons_sub_left_region_id_invalidC','select_neurons_sub_left_region_C','select_deact_neurons_sub_left_region_id_invalidC','select_deact_neurons_sub_left_region_C',...
    'select_neurons_sub_right_region_id_invalidC','select_neurons_sub_right_region_C','select_deact_neurons_sub_right_region_id_invalidC','select_deact_neurons_sub_right_region_C',...
    'deact_response_condition','-v7.3');

