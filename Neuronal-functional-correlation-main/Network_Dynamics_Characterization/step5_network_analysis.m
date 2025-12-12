clc; clear; close all; set(0,'defaultfigurecolor','w');
%% small world analysis for neurons network
% WLF 20230904

addpath(genpath('D:\PPT_gather\Neuronal functional connectivity\Utils'));

%% >>>>>>>>> load color
% colormap {light red; light blue; green; blue; light pink; light purple; light green; red; golden; light brown}
color_scheme_npg = get_color();
color_scheme_npg1 = color_scheme_npg(2:end,:); % the first color is red (don't use this)

%% ----------------------------------------------------------------------------------------------------------------------------------------------- show typical trace
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% load data by indicating experiment index
folder_path = ['D:\PPT_gather\Neuronal functional connectivity\dataset\'];
param.res_nums_thre = 8; % change here select neurons by response nums % adjust here 1
figure_path = [folder_path sprintf('result_for_thre%d\\',param.res_nums_thre)];
result_for_response_to_stis_path = [figure_path '\'];
mkdir(result_for_response_to_stis_path);

% region_data = [folder_path 'neuron_region_lr.mat'];
% load(region_data);

input_C_data = [result_for_response_to_stis_path 'neuron_region_lr_with_select_neurons.mat'];
load(input_C_data);
C = C_mat;

input_C_data = [result_for_response_to_stis_path 'neuron_region_lr_with_select_neurons_win12.mat'];
load(input_C_data);

network_analysis_for_select_neurons = 1;
network_analysis_for_all_neurons = 0;
network_analysis_for_sub_regions = 0;
network_analysis_for_select_sub_regions = 0;

param.data_session_n = 1;
param.edeg_threshold = 0.95

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% load stimuli
load([folder_path 'demo_visual_stimuli_with_label.mat']);
stimuli = load([folder_path 'demo_visual_stimuli_with_label.mat']);
stimuli_for_edge = load([folder_path 'demo_visual_stimuli_with_label.mat']);
stimuli_for_edge.start_edge = stimuli_for_edge.start_edge - 12;
stimuli_for_edge.end_edge = stimuli_for_edge.end_edge - 12;

% >>>>>> preprocess for different lable
buff = min(stimuli.stimuili_label_ind)-1;
stimuli.stimuili_label = stimuli.stimuili_label_ind'-buff;
stimuli.stimuili_label_ind = stimuli.stimuili_label_ind-buff;
stimuli.stimuli_array_with_label = stimuli.stimuli_array_with_label-buff;
stimuli.stimuli_array_with_label(stimuli.stimuli_array_with_label < 0) = 0;
tabulated_for_sti = tabulate(stimuli.stimuili_label_ind);
% >>>>>> process for null sti num
sti_num = size(tabulated_for_sti,1);
sti_num = find(tabulated_for_sti(:,2)~=0)
trial_nums = min(tabulated_for_sti(:,2));

buff_start = []; buff_end = [];
for i = 1:length(sti_num)
    trials_ind = find(stimuli.stimuili_label_ind==sti_num(i));
    buff_start = [buff_start;start_edge(trials_ind(1:trial_nums))];
    buff_end = [buff_end;end_edge(trials_ind(1:trial_nums))];
end
buff_start = buff_start'; buff_end =buff_end';
sti_num_vec = sti_num;
sti_num = length(sti_num);
% >>>>>>>> sti analysis only for train data
stimuli.stimuili_label = stimuli.stimuili_label(1:sti_num*trial_nums);
stimuli.stimuili_label_ind = stimuli.stimuili_label_ind(1:sti_num*trial_nums);
stimuli.start_edge = stimuli.start_edge(1:sti_num*trial_nums);
stimuli.end_edge = stimuli.end_edge(1:sti_num*trial_nums);

% ----------------------------------------------------------------------------------------------------------------------------------------------- network analysis （select neurons）
if network_analysis_for_select_neurons == 1

    % try
    [select_neurons_similarity_index,select_neurons_mean_similarity,...
        select_neurons_degree_sequence,select_neurons_closeness_sequence,select_neurons_betweenness_sequence,select_neurons_average_path_length_sequence,select_neurons_clustering_sequence,...
        select_neurons_all_degree,select_neurons_all_closeness,select_neurons_all_betweenness,select_neurons_all_average_path_length,select_neurons_all_clustering] = ...
    graph_network_analysis(C_mat,select_neurons_C_sort_by_region_edge_cell,stimuli,stimuli_for_edge,param);
    save([result_for_response_to_stis_path sprintf('index_statistic_all_trial_select_neurons%.2f.mat',param.edeg_threshold)],'param', ...
        'select_neurons_similarity_index','select_neurons_mean_similarity',...
        'select_neurons_degree_sequence','select_neurons_closeness_sequence','select_neurons_betweenness_sequence','select_neurons_average_path_length_sequence','select_neurons_clustering_sequence',...
        'select_neurons_all_degree','select_neurons_all_closeness','select_neurons_all_betweenness','select_neurons_all_average_path_length','select_neurons_all_clustering','-v7.3');
    %     catch
    %     end

    % % rich club analysis
    % for sti_i = 1:size(select_neurons_all_degree,1)
    %     for session_i = 1:size(select_neurons_all_degree,2)
    %         for trial_i = 1:size(select_neurons_all_degree,3)
    %             % for frame_i = 1:size(select_neurons_all_degree,4)
    %             for frame_i = 12
    %                 curr_brain_network_degree = select_neurons_all_degree(sti_i,session_i,trial_i,frame_i);
    %                 curr_brain_network_degree = curr_brain_network_degree{1};
    %                 rich_node_threshold = mean(curr_brain_network_degree);
    %                 richClubNodes_id = find(curr_brain_network_degree > rich_node_threshold);
    %             end
    %         end
    %     end
    % end

end


%% ----------------------------------------------------------------------------------------------------------------------------------------------- network analysis （all neurons）
if network_analysis_for_all_neurons == 1
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% load all neuron_corr
    %     load([result_for_response_to_stis_path 'all_neuron_corr.mat']);
    
    %     try
    [all_neurons_similarity_index,all_neurons_mean_similarity,...
        all_neurons_degree_sequence,all_neurons_closeness_sequence,all_neurons_betweenness_sequence,all_neurons_average_path_length_sequence,all_neurons_clustering_sequence,...
        all_neurons_all_degree,all_neurons_all_closeness,all_neurons_all_betweenness,all_neurons_all_average_path_length,all_neurons_all_clustering] = ...
    graph_network_analysis(C_mat,all_neurons_C_sort_by_region_edge_cell,stimuli,stimuli_for_edge,param);
    save([result_for_response_to_stis_path sprintf('index_statistic_all_trial_all_neurons_sti2_frame12%.2f.mat',param.edeg_threshold)],'param', ...
        'all_neurons_similarity_index','all_neurons_mean_similarity',...
        'all_neurons_degree_sequence','all_neurons_closeness_sequence','all_neurons_betweenness_sequence','all_neurons_average_path_length_sequence','all_neurons_clustering_sequence',...
        'all_neurons_all_degree','all_neurons_all_closeness','all_neurons_all_betweenness','all_neurons_all_average_path_length','all_neurons_all_clustering','-v7.3');
    %     catch
    %     end
end

%% ----------------------------------------------------------------------------------------------------------------------------------------------- network analysis (sub region)
if network_analysis_for_sub_regions == 1
    % try
    load([result_for_response_to_stis_path 'sub_region_corr.mat']);
    % get corr trials data for show
    sub_region_similarity_index = {};
    sub_region_mean_similarity = {};
    
    sub_region_degree_sequence = {};
    sub_region_closeness_sequence = {};
    sub_region_betweenness_sequence = {};
    sub_region_average_path_length_sequence = {};
    sub_region_clustering_sequence = {};
    
    sub_region_degree = {};
    sub_region_closeness = {};
    sub_region_betweenness = {};
    sub_region_average_path_length = {};
    sub_region_clustering = {};
    
    for regioni = 1:size(all_region_right_neurons_C_sort_by_region_edge_cell,2)
    %     try
        curr_region_right_neurons_C_sort_by_region_edge_cell = all_region_right_neurons_C_sort_by_region_edge_cell{1,regioni};
        [similarity_index,mean_similarity,degree_sequence,closeness_sequence,betweenness_sequence,average_path_length_sequence,clustering_sequence,...
            all_degree,all_closeness,all_betweenness,all_average_path_length,all_clustering] = ...
        graph_network_analysis(C_mat,curr_region_right_neurons_C_sort_by_region_edge_cell,stimuli,stimuli_for_edge,param);
        sub_region_similarity_index{regioni} = similarity_index;
        sub_region_mean_similarity{regioni} = mean_similarity;
    
        sub_region_degree_sequence{regioni} = degree_sequence;
        sub_region_closeness_sequence{regioni} = closeness_sequence;
        sub_region_betweenness_sequence{regioni} = betweenness_sequence;
        sub_region_average_path_length_sequence{regioni} = average_path_length_sequence;
        sub_region_clustering_sequence{regioni} = clustering_sequence;
    
        sub_region_degree{regioni} = degree_sequence;
        sub_region_closeness{regioni} = closeness_sequence;
        sub_region_betweenness{regioni} = betweenness_sequence;
        sub_region_average_path_length{regioni} = average_path_length_sequence;
        sub_region_clustering{regioni} = clustering_sequence;
    
    %     catch
    %     end
    end
    save([result_for_response_to_stis_path sprintf('index_statistic_all_trial_sub_region%.2f.mat',param.edeg_threshold)],'param', ...
        'sub_region_similarity_index','sub_region_mean_similarity',...
        'sub_region_degree_sequence','sub_region_closeness_sequence','sub_region_betweenness_sequence','sub_region_average_path_length_sequence','sub_region_clustering_sequence',...
        'sub_region_degree','sub_region_closeness','sub_region_betweenness','sub_region_average_path_length','sub_region_clustering','-v7.3');

end

%% ----------------------------------------------------------------------------------------------------------------------------------------------- network analysis (sub region)
if network_analysis_for_select_sub_regions == 1
    % try
    load([result_for_response_to_stis_path 'sub_region_corr_select.mat']);
    param.data_session_n = 1;
    % get corr trials data for show
    select_sub_region_similarity_index = {};
    select_sub_region_mean_similarity = {};
    
    select_sub_region_degree_sequence = {};
    select_sub_region_closeness_sequence = {};
    select_sub_region_betweenness_sequence = {};
    select_sub_region_average_path_length_sequence = {};
    select_sub_region_clustering_sequence = {};
    
    select_sub_region_degree = {};
    select_sub_region_closeness = {};
    select_sub_region_betweenness = {};
    select_sub_region_average_path_length = {};
    select_sub_region_clustering = {};
    
    for regioni = 1:size(all_region_select_right_neurons_C_sort_by_region_edge_cell,2)
        % try
        curr_all_region_select_right_neurons_C_sort_by_region_edge_celll = all_region_select_right_neurons_C_sort_by_region_edge_cell{1,regioni};
        [similarity_index,mean_similarity,degree_sequence,closeness_sequence,betweenness_sequence,average_path_length_sequence,clustering_sequence,...
            all_degree,all_closeness,all_betweenness,all_average_path_length,all_clustering] = ...
        graph_network_analysis(C_mat,all_region_select_right_neurons_C_sort_by_region_edge_cell{regioni},stimuli,stimuli_for_edge,param);
        select_sub_region_similarity_index{regioni} = similarity_index;
        select_sub_region_mean_similarity{regioni} = mean_similarity;
    
        select_sub_region_degree_sequence{regioni} = degree_sequence;
        select_sub_region_closeness_sequence{regioni} = closeness_sequence;
        select_sub_region_betweenness_sequence{regioni} = betweenness_sequence;
        select_sub_region_average_path_length_sequence{regioni} = average_path_length_sequence;
        select_sub_region_clustering_sequence{regioni} = clustering_sequence;
    
        select_sub_region_degree{regioni} = degree_sequence;
        select_sub_region_closeness{regioni} = closeness_sequence;
        select_sub_region_betweenness{regioni} = betweenness_sequence;
        select_sub_region_average_path_length{regioni} = average_path_length_sequence;
        select_sub_region_clustering{regioni} = clustering_sequence;
    
        % catch
        % end
    end

    save([result_for_response_to_stis_path sprintf('index_statistic_all_trial_select_sub_region%.2f.mat',param.edeg_threshold)],'param', ...
        'select_sub_region_similarity_index','select_sub_region_mean_similarity',...
        'select_sub_region_degree_sequence','select_sub_region_closeness_sequence',...
        'select_sub_region_betweenness_sequence','select_sub_region_average_path_length_sequence',...
        'select_sub_region_clustering_sequence','select_sub_region_degree',...
        'select_sub_region_closeness','select_sub_region_betweenness'...
        ,'select_sub_region_average_path_length','select_sub_region_clustering','-v7.3');

end


