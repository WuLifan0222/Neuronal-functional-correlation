clc; clear; close all; set(0,'defaultfigurecolor','w');
%% quantify trial-to-trial variability
% WLF 20241220
addpath(genpath('D:\PPT_gather\Neuronal functional connectivity\Utils'));

color_scheme_npg = get_color();
red_blue_color = get_red_blue_color_map();

%% >>> choose data sti kind
data = '' 
% data = 'drifting' 
% data = 'natureimage' 

%% >>> set param
top_n = 200;
process_nozeros = 1;

%% >>>>>>>>> Load experiment (load data and stimuli)
folder_path = ['D:\PPT_gather\Neuronal functional connectivity\dataset\'];
param.res_nums_thre = 8; % change here select neurons by response nums % adjust here 1
figure_path = [folder_path sprintf('result_for_thre%d\\',param.res_nums_thre)];
result_for_response_to_stis_path = [figure_path '\'];
mkdir(result_for_response_to_stis_path);

% load  all edge
load([figure_path 'neuron_region_lr_with_select_neurons.mat']);
stimuli = load([folder_path 'demo_visual_stimuli_with_label.mat']);
load([figure_path 'neuron_region_lr_with_select_neurons.mat']);
% load([figure_path 'top_100_edges_trial_variance.mat']);

%% >>>>>>>>> atlas
atlas_bg = importdata('D:\PPT_gather\Neuronal functional connectivity\Utils\atlas\cortical_out_line_resize_5_30.mat');
atlas = importdata('D:\PPT_gather\Neuronal functional connectivity\Utils\atlas\atlas_top_projection_30.mat');
resize_ratio = 5;

%% get trial variability value for all edge & neurons
% >>>>>>>>>>>>>>>>>>>>>>> TV select from all neuron
[all_neuron_stdResult,sort_all_neuron_stdResult_id,sort_all_neuron_stdResult_value] = get_all_stdResult_value(C_mat,stimuli, param);
compare_sort_all_neuron_stdResult_value{1} = sort_all_neuron_stdResult_value{5}; %blank
sti_sort_all_neuron_stdResult_value0 = [sort_all_neuron_stdResult_value{2},...
                                            sort_all_neuron_stdResult_value{3},...
                                            sort_all_neuron_stdResult_value{4},...
                                            sort_all_neuron_stdResult_value{1}];
compare_sort_all_neuron_stdResult_value{2} = mean(sti_sort_all_neuron_stdResult_value0,2); % sti average


% % >>>>>>>>>>>>>>>>>>>>>>> ttest neurons
[select_neuron_stdResult,sort_select_neuron_stdResult_id,sort_select_neuron_stdResult_value] = get_all_stdResult_value(select_neurons_C,stimuli, param);
compare_sort_select_neuron_stdResult_value{1} = sort_select_neuron_stdResult_value{5}; %blank
sti_sort_select_neuron_stdResult_value0 = [sort_select_neuron_stdResult_value{2},...
                                            sort_select_neuron_stdResult_value{3},...
                                            sort_select_neuron_stdResult_value{4},...
                                            sort_select_neuron_stdResult_value{1}];
compare_sort_select_neuron_stdResult_value{2} = mean(sti_sort_select_neuron_stdResult_value0,2); % sti average


% % >>>>>>>>>>>>>>>>>>>>>>> TV select edge
[edge_stdResult,sort_edge_stdResult_id,sort_edge_stdResult_value] = get_all_stdResult_value(select_neurons_C_sort_by_region_corr_matrix,stimuli, param);
compare_sort_edge_stdResult_value{1} = sort_edge_stdResult_value{5}; %blank
sti_sort_edge_stdResult_value0 = [sort_edge_stdResult_value{2},...
                                            sort_edge_stdResult_value{3},...
                                            sort_edge_stdResult_value{4},...
                                            sort_edge_stdResult_value{1}];
compare_sort_edge_stdResult_value{2} = mean(sti_sort_edge_stdResult_value0,2); % sti average
select_neuron_num = size(select_neurons_C,1);


if strcmp(data, 'drifting')
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% edge
    for i = 1:10
        top_values = [];
        top_trace_id = [];
        for sti_i = 1:4
            % 求每组长度
            group_size = floor(length(sort_edge_stdResult_value{sti_i}) / 10);
            % 每组的起始和结束索引
            start_idx = (i-1)*group_size + 1;
            if i < 10
                end_idx = i*group_size;
            else
                % 最后一组包含剩余所有元素
                end_idx = length(sort_edge_stdResult_value{sti_i});
            end
            % 取当前组的数据
            group_data_value = sort_edge_stdResult_value{sti_i}(start_idx:end_idx);
            group_data_value_id = sort_edge_stdResult_id{sti_i}(start_idx:end_idx);

            % 取前n个最大值
            curr_sti_top_values = group_data_value(1:min(top_n, length(sort_edge_stdResult_value{sti_i})));
            curr_sti_top_trace_id = group_data_value_id(1:min(top_n, length(sort_edge_stdResult_value{sti_i})));
            top_values = [top_values;curr_sti_top_values];
            top_trace_id = [top_trace_id;curr_sti_top_trace_id];
        end
        top_values_mean = mean(top_values,1);
        curr_stage_edges = select_neurons_C_sort_by_region_corr_matrix(top_trace_id,:);
        deconve_trace_with_label_gray(curr_stage_edges, stimuli, color_scheme_npg);
        savefig([figure_path sprintf('stage%d_edges_show_deconve',i)]);
        exportgraphics(gcf,[figure_path sprintf('stage%d_edges_show_deconve',i),'.png'],'Resolution',300)
        save([figure_path 'edge_stage' num2str(i) '.mat'], 'top_values','top_values_mean','top_trace_id','curr_stage_edges','-v7.3');
    end

    % add shuffle data compare
    top_values = [];
    top_trace_id = [];
    for sti_i = 1:4
        random_numbers = randi([1, length(sort_edge_stdResult_value{sti_i})], top_n, 1);  
        curr_sti_top_values = sort_edge_stdResult_value{sti_i}(random_numbers);
        curr_sti_top_trace_id = sort_edge_stdResult_id{sti_i}(random_numbers);
        top_values = [top_values;curr_sti_top_values];
        top_trace_id = [top_trace_id;curr_sti_top_trace_id];
    end
    top_values_mean = mean(top_values,1);
    curr_stage_edges = select_neurons_C_sort_by_region_corr_matrix(top_trace_id,:);
    deconve_trace_with_label_gray(curr_stage_edges, stimuli, color_scheme_npg);
    savefig([figure_path sprintf('shuffle_edges_show_deconve',i)]);
    exportgraphics(gcf,[figure_path sprintf('shuffle_edges_show_deconve',i),'.png'],'Resolution',300)
    save([figure_path 'edge_shuffle'  '.mat'], 'top_values','top_values_mean','top_trace_id','curr_stage_edges','-v7.3');

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % neuron
    for i = 1:10
    top_values = [];
    top_trace_id = [];
    for sti_i = 1:4
        % 求每组长度
        group_size = floor(length(sort_all_neuron_stdResult_value{sti_i}) / 10);
        % 每组的起始和结束索引
        start_idx = (i-1)*group_size + 1;
        if i < 10
            end_idx = i*group_size;
        else
            % 最后一组包含剩余所有元素
            end_idx = length(sort_all_neuron_stdResult_value{sti_i});
        end
        % 取当前组的数据
        group_data_value = sort_all_neuron_stdResult_value{sti_i}(start_idx:end_idx);
        group_data_value_id = sort_all_neuron_stdResult_id{sti_i}(start_idx:end_idx);

        % 取前top_n个最大值
        curr_sti_top_values = group_data_value(1:min(top_n, length(sort_all_neuron_stdResult_value{sti_i})));
        curr_sti_top_trace_id = group_data_value_id(1:min(top_n, length(sort_all_neuron_stdResult_value{sti_i})));
        top_values = [top_values;curr_sti_top_values];
        top_trace_id = [top_trace_id;curr_sti_top_trace_id];
    end
    top_values_mean = mean(top_values,1);
    curr_stage_all_neurons = C_mat(top_trace_id,:);
    deconve_trace_with_label_gray(curr_stage_all_neurons, stimuli, color_scheme_npg);
    savefig([figure_path sprintf('stage%d_all_neurons_show_deconve',i)]);
    exportgraphics(gcf,[figure_path sprintf('stage%d_all_neurons_show_deconve',i),'.png'],'Resolution',300)
    save([figure_path 'all_neuron_stage' num2str(i) '.mat'], 'top_values','top_values_mean','top_trace_id','curr_stage_all_neurons','-v7.3');
end

% add shuffle data compare
top_values = [];
top_trace_id = [];
for sti_i = 1:4
    random_numbers = randi([1, length(sort_all_neuron_stdResult_value{sti_i})], top_n, 1);  
    curr_sti_top_values = sort_all_neuron_stdResult_value{sti_i}(random_numbers);
    curr_sti_top_trace_id = sort_all_neuron_stdResult_id{sti_i}(random_numbers);
    top_values = [top_values;curr_sti_top_values];
    top_trace_id = [top_trace_id;curr_sti_top_trace_id];
end
top_values_mean = mean(top_values,1);
curr_stage_all_neurons = C_mat(top_trace_id,:);
deconve_trace_with_label_gray(curr_stage_all_neurons, stimuli, color_scheme_npg);
savefig([figure_path sprintf('shuffle_all_neurons_show_deconve',i)]);
exportgraphics(gcf,[figure_path sprintf('shuffle_all_neurons_show_deconve',i),'.png'],'Resolution',300)
save([figure_path 'all_neuron_shuffle'  '.mat'], 'top_values','top_values_mean','top_trace_id','curr_stage_all_neurons','-v7.3');

end

if strcmp(data, 'natureimage')
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% edge
    for i = 1:10
        top_values = [];
        top_trace_id = [];
        for sti_i = 2:6
            % 求每组长度
            group_size = floor(length(sort_edge_stdResult_value{sti_i}) / 10);
            % 每组的起始和结束索引
            start_idx = (i-1)*group_size + 1;
            if i < 10
                end_idx = i*group_size;
            else
                % 最后一组包含剩余所有元素
                end_idx = length(sort_edge_stdResult_value{sti_i});
            end
            % 取当前组的数据
            group_data_value = sort_edge_stdResult_value{sti_i}(start_idx:end_idx);
            group_data_value_id = sort_edge_stdResult_id{sti_i}(start_idx:end_idx);

            % 取前top_n个最大值
            curr_sti_top_values = group_data_value(1:min(top_n, length(sort_edge_stdResult_value{sti_i})));
            curr_sti_top_trace_id = group_data_value_id(1:min(top_n, length(sort_edge_stdResult_value{sti_i})));
            top_values = [top_values;curr_sti_top_values];
            top_trace_id = [top_trace_id;curr_sti_top_trace_id];
        end
        top_values_mean = mean(top_values,1);
        curr_stage_edges = select_neurons_C_sort_by_region_corr_matrix(top_trace_id,:);
        deconve_trace_with_label_gray(curr_stage_edges, stimuli, color_scheme_npg);
        savefig([figure_path sprintf('stage%d_edges_show_deconve',i)]);
        exportgraphics(gcf,[figure_path sprintf('stage%d_edges_show_deconve',i),'.png'],'Resolution',300)
        save([figure_path 'edge_stage' num2str(i) '.mat'], 'top_values','top_values_mean','top_trace_id','curr_stage_edges','-v7.3');
    end

    % add shuffle data compare
    top_values = [];
    top_trace_id = [];
    for sti_i = 2:6
        random_numbers = randi([1, length(sort_edge_stdResult_value{sti_i})], top_n, 1);  
        curr_sti_top_values = sort_edge_stdResult_value{sti_i}(random_numbers);
        curr_sti_top_trace_id = sort_edge_stdResult_id{sti_i}(random_numbers);
        top_values = [top_values;curr_sti_top_values];
        top_trace_id = [top_trace_id;curr_sti_top_trace_id];
    end
    top_values_mean = mean(top_values,1);
    curr_stage_edges = select_neurons_C_sort_by_region_corr_matrix(top_trace_id,:);
    deconve_trace_with_label_gray(curr_stage_edges, stimuli, color_scheme_npg);
    savefig([figure_path sprintf('shuffle_edges_show_deconve',i)]);
    exportgraphics(gcf,[figure_path sprintf('shuffle_edges_show_deconve',i),'.png'],'Resolution',300)
    save([figure_path 'edge_shuffle'  '.mat'], 'top_values','top_values_mean','top_trace_id','curr_stage_edges','-v7.3');

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % neuron
    for i = 1:10
    top_values = [];
    top_trace_id = [];
    for sti_i = 2:6
        % 求每组长度
        group_size = floor(length(sort_all_neuron_stdResult_value{sti_i}) / 10);
        % 每组的起始和结束索引
        start_idx = (i-1)*group_size + 1;
        if i < 10
            end_idx = i*group_size;
        else
            % 最后一组包含剩余所有元素
            end_idx = length(sort_all_neuron_stdResult_value{sti_i});
        end
        % 取当前组的数据
        group_data_value = sort_all_neuron_stdResult_value{sti_i}(start_idx:end_idx);
        group_data_value_id = sort_all_neuron_stdResult_id{sti_i}(start_idx:end_idx);

        % 取前top_n个最大值
        curr_sti_top_values = group_data_value(1:min(top_n, length(sort_all_neuron_stdResult_value{sti_i})));
        curr_sti_top_trace_id = group_data_value_id(1:min(top_n, length(sort_all_neuron_stdResult_value{sti_i})));
        top_values = [top_values;curr_sti_top_values];
        top_trace_id = [top_trace_id;curr_sti_top_trace_id];
    end
    top_values_mean = mean(top_values,1);
    curr_stage_all_neurons = C_mat(top_trace_id,:);
    deconve_trace_with_label_gray(curr_stage_all_neurons, stimuli, color_scheme_npg);
    savefig([figure_path sprintf('stage%d_all_neurons_show_deconve',i)]);
    exportgraphics(gcf,[figure_path sprintf('stage%d_all_neurons_show_deconve',i),'.png'],'Resolution',300)
    save([figure_path 'all_neuron_stage' num2str(i) '.mat'], 'top_values','top_values_mean','top_trace_id','curr_stage_all_neurons','-v7.3');
end

% add shuffle data compare
top_values = [];
top_trace_id = [];
for sti_i = 2:6
    random_numbers = randi([1, length(sort_all_neuron_stdResult_value{sti_i})], top_n, 1);  
    curr_sti_top_values = sort_all_neuron_stdResult_value{sti_i}(random_numbers);
    curr_sti_top_trace_id = sort_all_neuron_stdResult_id{sti_i}(random_numbers);
    top_values = [top_values;curr_sti_top_values];
    top_trace_id = [top_trace_id;curr_sti_top_trace_id];
end
top_values_mean = mean(top_values,1);
curr_stage_all_neurons = C_mat(top_trace_id,:);
deconve_trace_with_label_gray(curr_stage_all_neurons, stimuli, color_scheme_npg);
savefig([figure_path sprintf('shuffle_all_neurons_show_deconve',i)]);
exportgraphics(gcf,[figure_path sprintf('shuffle_all_neurons_show_deconve',i),'.png'],'Resolution',300)
save([figure_path 'all_neuron_shuffle'  '.mat'], 'top_values','top_values_mean','top_trace_id','curr_stage_all_neurons','-v7.3');

end



%% show select edge in atlas

data_session_n = 1;
each_session_trials_num = 20;

% show musk
% % >>>>>> auto get more param
% sort for better visual
[~,~,select_neurons_C_sort_by_region_corr_matrix_visualsort,~, ~, ~, ~,~, ~] = sort_multi_neurons_trials_matrix_by_sti(select_neurons_C_sort_by_region_corr_matrix, stimuli, param);
% get corr trials data for show
[corr_trials_cell,before,during,after,sti_num] = get_corr_trials_matrix_by_sti(select_neurons_C_sort_by_region_edge_cell,stimuli,param);
C_matNorm = (C_mat - min(C_mat(:))) / (max(C_mat(:)) - min(C_mat(:)));
[C_matNorm_cell,before,during,after,sti_num] = get_trials_matrix_multi_trace_by_sti((C_matNorm+0.0000001),stimuli,param);
select_neurons_CNorm = (select_neurons_C - min(select_neurons_C(:))) / (max(select_neurons_C(:)) - min(select_neurons_C(:)));
[select_neurons_CNorm_cell,before,during,after,sti_num] = get_trials_matrix_multi_trace_by_sti((select_neurons_CNorm+0.0000001),stimuli,param);


for sti_ii = 1:param.sti_num
    curr_select_edge_id = sort_edge_stdResult_id{sti_ii}(1:select_neuron_num); 
    curr_select_edge_visualsort = select_neurons_C_sort_by_region_corr_matrix_visualsort(curr_select_edge_id,:);
    each_sti_edge_select_by_trials_sort{sti_ii} = select_neurons_C_sort_by_region_corr_matrix(curr_select_edge_id,:);
    point_num = size(select_neurons_C,1);
    [curr_select_sti_musk] = select_edge_transfer_to_musk(curr_select_edge_id,select_neurons_C_sort_by_region_corr_matrix_visualsort,point_num);
    select_sti_musk{sti_ii} = curr_select_sti_musk;
end

mean_R_cell = {};
data_sessioni = data_session_n;
for sti_i = 1:param.sti_num
    curr_sti_R = corr_trials_cell{sti_i};
    curr_sti_all_trials_R = [];
    mean_sti_R = [];
    for trial_frame_i = 1:size(curr_sti_R,2)
        mean_sti_frame_R = [];
        for trial_i = each_session_trials_num*(data_sessioni-1)+1:each_session_trials_num*data_sessioni
            curr_sti_trial_frame_R = curr_sti_R{trial_i,trial_frame_i};
            curr_sti_frame_R(trial_i,:,:) = curr_sti_trial_frame_R;
        end
        point_num = size(curr_sti_trial_frame_R,1);
        mean_sti_frame_R = mean(curr_sti_frame_R,1);
        mean_sti_frame_R = reshape(mean_sti_frame_R,[point_num,point_num]);
        mean_sti_R{trial_frame_i} = mean_sti_frame_R;
    end
    mean_R_cell{sti_i} = mean_sti_R;
end

point_num = size(select_neurons_C,1);
for sti_ii = 1:param.sti_num
    curr_corr_sti = mean_R_cell{sti_ii};% curr_corr_sti = mean_R_cell{sti_iii};

    for frame_i = 1:(before+during+after)
    curr_corr_sti_trial{frame_i} = curr_corr_sti{1,frame_i};
    end

    % 1 means show first trial, 12 is key frame
    R = curr_corr_sti_trial{12};
    edge_thre= 0.5;
    C_value= select_neurons_CNorm_cell{sti_ii}{1}(:,12);% C_value= select_neurons_CNorm_cell{sti_iii}{trial_i}(:,corr_i) ;
    C_all_value = C_matNorm_cell{sti_ii}{1}(:,12) ;% C_all_value = C_matNorm_cell{sti_iii}{trial_i}(:,corr_i) ;
    text_title = sprintf('sti%d trial 1 during sti',sti_ii);
    curr_select_sti_musk = select_sti_musk{sti_ii};% curr_select_sti_musk = select_sti_musk{sti_iii};
    F_musk = R.* curr_select_sti_musk;
    % figure;imagesc(F_musk)
    try
    figure
    plot_single_frame_cluster_connection_in_atlas(C_mat,F_musk,F_musk,edge_thre,C_all_value,C_value,select_neurons_id,atlas,xlabel_cell,brain_region,neuron_position_x_mat,neuron_position_y_mat,color_scheme_npg(sti_ii+1,:),text_title)
    hold on; axis([800, 1; 2100, 1200]/5); 
    axis off
    savefig([figure_path sprintf('TV_select_top_%d_stable_edge_position_sti_%d',select_neuron_num,sti_ii)]);
    exportgraphics(gcf,[figure_path sprintf('TV_select_top_%d_stable_edge_position_sti_%d',select_neuron_num,sti_ii),'.png'],'Resolution',300)
    catch
    end
end
close all

for show_std_sti_i = 1:size(sort_all_neuron_stdResult_id,2)

%% statistic TV for sub regions

% define edge kind
TV_value_within_VISp = [];
TV_value_VISp_other = [];
TV_value_within_other = [];
select_neuron_region_id = neuron_region_id(select_neurons_id);

targetString = 'VISp1'; 
isMatch = strcmp(xlabel_annoy, targetString);
VISp_id = find(isMatch == 1);

curr_select_edge_id = sort_edge_stdResult_id{show_std_sti_i}(1:select_neuron_num); 
curr_select_edge_TV_value = sort_edge_stdResult_value{show_std_sti_i}(1:select_neuron_num); 
for select_top_edge_i =1:size(curr_select_edge_id,1)
    curr_id = curr_select_edge_id(select_top_edge_i);
    curr_TV_value = curr_select_edge_TV_value(select_top_edge_i);
    [curr_musk] = select_edge_transfer_to_musk(curr_id,select_neurons_C_sort_by_region_corr_matrix_visualsort,point_num);
    [curr_edge_i,curr_edge_j] = find(curr_musk~=0)
    crrr_edge_region_i = select_neuron_region_id(curr_edge_i) ;
    crrr_edge_region_j = select_neuron_region_id(curr_edge_j) ;
    %  >>>>>>>> within VISp
    if ((crrr_edge_region_i==VISp_id)&&(crrr_edge_region_j==VISp_id))    
        TV_value_within_VISp = [TV_value_within_VISp,curr_TV_value];
    %  >>>>>>>>  VISp & other region
    elseif ((crrr_edge_region_i==VISp_id)||(crrr_edge_region_j==VISp_id))
        TV_value_VISp_other = [TV_value_VISp_other,curr_TV_value];
    %  >>>>>>>> within other region
    elseif ((crrr_edge_region_i~=VISp_id)&&(crrr_edge_region_j~=VISp_id))
        TV_value_within_other = [TV_value_within_other,curr_TV_value];
    end
end


% define neuron kind
TV_neuron_value_within_VISp = [];
TV_neuron_value_VISp_other = [];
TV_neuron_value_within_other = [];

select_neuron_region_id = neuron_region_id(select_neurons_id);
targetString = 'VISp1'; 
isMatch = strcmp(xlabel_annoy, targetString);
VISp_id = find(isMatch == 1);

curr_select_neuron_id = sort_select_neuron_stdResult_id{1}(1:select_neuron_num); 
curr_select_neuron_TV_value = sort_select_neuron_stdResult_value{1}(1:select_neuron_num); 
for select_top_neuron_i =1:size(curr_select_neuron_id,1)
    curr_neuron_id = curr_select_neuron_id(select_top_neuron_i);
    curr_TV_value = curr_select_neuron_TV_value(select_top_neuron_i);
    curr_neuron_region_i = select_neuron_region_id(select_top_neuron_i) ;
    %  >>>>>>>> within VISp
    if ((curr_neuron_region_i==VISp_id))    
        TV_neuron_value_within_VISp = [TV_neuron_value_within_VISp,curr_TV_value];
    %  >>>>>>>> within other region
    elseif ((curr_neuron_region_i~=VISp_id))
        TV_neuron_value_within_other = [TV_neuron_value_within_other,curr_TV_value];
    end
end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% VISp VIS SS MO RSP
%% statistic TV for sub regions

% define edge kind
TV_value_within_VISp = [];
TV_value_within_VIS = [];
TV_value_within_SS = [];
TV_value_within_RSP = [];
TV_value_within_MO = [];
TV_value_within_AUD = [];

TV_value_within_VISa = [];
TV_value_within_VISal = [];
TV_value_within_VISam = [];
TV_value_within_VISl = [];
TV_value_within_VISli = [];
TV_value_within_VISpl = [];
TV_value_within_VISpor = [];
TV_value_within_VISrl = [];
TV_value_within_VISpm = [];

select_neuron_region_id = neuron_region_id(select_neurons_id);

targetString = 'VISp1'; 
isMatch = strcmp(xlabel_annoy, targetString);
VISp_id = find(isMatch == 1);

target_labels = {'VISa1','VISal1','VISam1','VISl1','VISli1','VISpl1','VISpor1','VISrl1','VISpm1'};
isMatch = ismember(xlabel_annoy, target_labels);
VIS_id = find(isMatch == 1);

target_labels = {'VISa1'};
isMatch = ismember(xlabel_annoy, target_labels);
VISa_id = find(isMatch == 1);

target_labels = {'VISal1'};
isMatch = ismember(xlabel_annoy, target_labels);
VISal_id = find(isMatch == 1);

target_labels = {'VISam1'};
isMatch = ismember(xlabel_annoy, target_labels);
VISam_id = find(isMatch == 1);

target_labels = {'VISl1'};
isMatch = ismember(xlabel_annoy, target_labels);
VISl_id = find(isMatch == 1);

target_labels = {'VISli1'};
isMatch = ismember(xlabel_annoy, target_labels);
VISli_id = find(isMatch == 1);

target_labels = {'VISpl1'};
isMatch = ismember(xlabel_annoy, target_labels);
VISpl_id = find(isMatch == 1);

target_labels = {'VISpor1'};
isMatch = ismember(xlabel_annoy, target_labels);
VISpor_id = find(isMatch == 1);

target_labels = {'VISrl1'};
isMatch = ismember(xlabel_annoy, target_labels);
VISrl_id = find(isMatch == 1);

target_labels = {'VISpm1'};
isMatch = ismember(xlabel_annoy, target_labels);
VISpm_id = find(isMatch == 1);

target_labels = {'SSp-bfd1','SSp-ll1','SSp-m1','SSp-n1','SSp-tr1','SSp-un1','SSp-ul1','SSs1'};
isMatch = ismember(xlabel_annoy, target_labels);
SS_id = find(isMatch == 1);

target_labels = {'RSPd1','RSPagl1'};
isMatch = ismember(xlabel_annoy, target_labels);
RSP_id = find(isMatch == 1);

target_labels = {'MOp1'};
isMatch = ismember(xlabel_annoy, target_labels);
MO_id = find(isMatch == 1);

target_labels = {'AUDd1','AUDpo1','AUDp1'};
isMatch = ismember(xlabel_annoy, target_labels);
AUD_id = find(isMatch == 1);

curr_select_edge_id = sort_edge_stdResult_id{show_std_sti_i}(1:select_neuron_num); 
curr_select_edge_TV_value = sort_edge_stdResult_value{show_std_sti_i}(1:select_neuron_num); 

for select_top_edge_i =1:size(curr_select_edge_id,1)
    curr_id = curr_select_edge_id(select_top_edge_i);
    curr_TV_value = curr_select_edge_TV_value(select_top_edge_i);
    [curr_musk] = select_edge_transfer_to_musk(curr_id,select_neurons_C_sort_by_region_corr_matrix_visualsort,point_num);
    [curr_edge_i,curr_edge_j] = find(curr_musk~=0)
    crrr_edge_region_i = select_neuron_region_id(curr_edge_i) ;
    crrr_edge_region_j = select_neuron_region_id(curr_edge_j) ;
    %  >>>>>>>> within VISp
    if ((crrr_edge_region_i==VISp_id)||(crrr_edge_region_j==VISp_id))    
        TV_value_within_VISp = [TV_value_within_VISp,curr_TV_value];   
    %  >>>>>>>>  VIS & other region
    elseif (ismember(crrr_edge_region_i,VIS_id)||ismember(crrr_edge_region_j,VIS_id))
        TV_value_within_VIS = [TV_value_within_VIS,curr_TV_value];
    elseif (ismember(crrr_edge_region_i,SS_id)||ismember(crrr_edge_region_j,SS_id))
        TV_value_within_SS = [TV_value_within_SS,curr_TV_value];
    elseif (ismember(crrr_edge_region_i,RSP_id)||ismember(crrr_edge_region_j,RSP_id))
        TV_value_within_RSP = [TV_value_within_RSP,curr_TV_value];
    elseif (ismember(crrr_edge_region_i,MO_id)||ismember(crrr_edge_region_j,MO_id))
        TV_value_within_MO = [TV_value_within_MO,curr_TV_value];
    elseif (ismember(crrr_edge_region_i,AUD_id)||ismember(crrr_edge_region_j,AUD_id))
        TV_value_within_AUD = [TV_value_within_AUD,curr_TV_value];
    end
end

for select_top_edge_i =1:size(curr_select_edge_id,1)
    curr_id = curr_select_edge_id(select_top_edge_i);
    curr_TV_value = curr_select_edge_TV_value(select_top_edge_i);
    [curr_musk] = select_edge_transfer_to_musk(curr_id,select_neurons_C_sort_by_region_corr_matrix_visualsort,point_num);
    [curr_edge_i,curr_edge_j] = find(curr_musk~=0)

    crrr_edge_region_i = select_neuron_region_id(curr_edge_i) ;
    crrr_edge_region_j = select_neuron_region_id(curr_edge_j) ;

    %  >>>>>>>> within VIS
    if (ismember(crrr_edge_region_i,VISa_id)||ismember(crrr_edge_region_j,VISa_id))
        TV_value_within_VISa = [TV_value_within_VISa,curr_TV_value];
    elseif (ismember(crrr_edge_region_i,VISal_id)||ismember(crrr_edge_region_j,VISal_id))
        TV_value_within_VISal = [TV_value_within_VISal,curr_TV_value];
    elseif (ismember(crrr_edge_region_i,VISam_id)||ismember(crrr_edge_region_j,VISam_id))
        TV_value_within_VISam = [TV_value_within_VISam,curr_TV_value];
    elseif (ismember(crrr_edge_region_i,VISl_id)||ismember(crrr_edge_region_j,VISl_id))
        TV_value_within_VISl = [TV_value_within_VISl,curr_TV_value];
    elseif (ismember(crrr_edge_region_i,VISli_id)||ismember(crrr_edge_region_j,VISli_id))
        TV_value_within_VISli = [TV_value_within_VISli,curr_TV_value];
    elseif (ismember(crrr_edge_region_i,VISpl_id)||ismember(crrr_edge_region_j,VISpl_id))
        TV_value_within_VISpl = [TV_value_within_VISpl,curr_TV_value];
    elseif (ismember(crrr_edge_region_i,VISpor_id)||ismember(crrr_edge_region_j,VISpor_id))
        TV_value_within_VISpor = [TV_value_within_VISpor,curr_TV_value];
    elseif (ismember(crrr_edge_region_i,VISrl_id)||ismember(crrr_edge_region_j,VISrl_id))
        TV_value_within_VISrl = [TV_value_within_VISrl,curr_TV_value];
    elseif (ismember(crrr_edge_region_i,VISpm_id)||ismember(crrr_edge_region_j,VISpm_id))
        TV_value_within_VISpm = [TV_value_within_VISpm,curr_TV_value];
    end

end


% define neuron kind
TV_neuron_value_within_VISp = [];
TV_neuron_value_within_VIS = [];
TV_neuron_value_within_SS = [];
TV_neuron_value_within_RSP = [];
TV_neuron_value_within_MO = [];
TV_neuron_value_within_AUD = [];

TV_neuron_value_within_VISa = [];
TV_neuron_value_within_VISal = [];
TV_neuron_value_within_VISam = [];
TV_neuron_value_within_VISl = [];
TV_neuron_value_within_VISli = [];
TV_neuron_value_within_VISpl = [];
TV_neuron_value_within_VISpor = [];
TV_neuron_value_within_VISrl = [];
TV_neuron_value_within_VISpm = [];

curr_select_neuron_id = sort_select_neuron_stdResult_id{1}(1:select_neuron_num); 
curr_select_neuron_TV_value = sort_select_neuron_stdResult_value{1}(1:select_neuron_num); 
for select_top_neuron_i =1:size(curr_select_neuron_id,1)
    curr_neuron_id = curr_select_neuron_id(select_top_neuron_i);
    curr_TV_value = curr_select_neuron_TV_value(select_top_neuron_i);
    curr_neuron_region_i = select_neuron_region_id(select_top_neuron_i) ;
    %  >>>>>>>> within VISp
    if ((curr_neuron_region_i==VISp_id))   
        TV_neuron_value_within_VISp = [TV_neuron_value_within_VISp,curr_TV_value];
    %  >>>>>>>>  VIS & other region
    elseif (ismember(curr_neuron_region_i,VIS_id))
        TV_neuron_value_within_VIS = [TV_neuron_value_within_VIS,curr_TV_value];
    elseif (ismember(curr_neuron_region_i,SS_id))
        TV_neuron_value_within_SS = [TV_neuron_value_within_SS,curr_TV_value];
    elseif (ismember(curr_neuron_region_i,RSP_id))
        TV_neuron_value_within_RSP = [TV_neuron_value_within_RSP,curr_TV_value];
    elseif (ismember(curr_neuron_region_i,MO_id))
        TV_neuron_value_within_MO = [TV_neuron_value_within_MO,curr_TV_value];
    elseif (ismember(curr_neuron_region_i,AUD_id))
        TV_neuron_value_within_AUD = [TV_value_within_AUD,curr_TV_value];
    end
end

for select_top_neuron_i =1:size(curr_select_neuron_id,1)
    curr_neuron_id = curr_select_neuron_id(select_top_neuron_i);
    curr_TV_value = curr_select_neuron_TV_value(select_top_neuron_i);
    curr_neuron_region_i = select_neuron_region_id(select_top_neuron_i) ;
    %  >>>>>>>> within VIS
    if ((curr_neuron_region_i==VISa_id))   
        TV_neuron_value_within_VISa = [TV_neuron_value_within_VISa,curr_TV_value];
    elseif (ismember(curr_neuron_region_i,VISal_id))
        TV_neuron_value_within_VISal = [TV_neuron_value_within_VISal,curr_TV_value];
    elseif (ismember(curr_neuron_region_i,VISam_id))
        TV_neuron_value_within_VISam = [TV_neuron_value_within_VISam,curr_TV_value];
    elseif (ismember(curr_neuron_region_i,VISl_id))
        TV_neuron_value_within_VISl = [TV_neuron_value_within_VISl,curr_TV_value];
    elseif (ismember(curr_neuron_region_i,VISli_id))
        TV_neuron_value_within_VISli = [TV_neuron_value_within_VISli,curr_TV_value];
    elseif (ismember(curr_neuron_region_i,VISpl_id))
        TV_neuron_value_within_VISpl = [TV_neuron_value_within_VISpl,curr_TV_value];
    elseif (ismember(curr_neuron_region_i,VISpor_id))
        TV_neuron_value_within_VISpor = [TV_neuron_value_within_VISpor,curr_TV_value];
    elseif (ismember(curr_neuron_region_i,VISrl_id))
        TV_neuron_value_within_VISrl = [TV_neuron_value_within_VISrl,curr_TV_value];
    elseif (ismember(curr_neuron_region_i,VISpm_id))
        TV_neuron_value_within_VISpm = [TV_neuron_value_within_VISpm,curr_TV_value];
    end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% jitter plot with significance test

% 定义脑区列表
region_names = {'VISp', 'VIS', 'VISa','VISal','VISam','VISl','VISli','VISpl','VISpor','VISrl','VISpm','SS', 'RSP', 'MO', 'AUD'};

% 抖动参数
jitter_amount = 0.15; 
neuron_color = [0.2 0.5 0.2];
edge_color = [1 0.85 0.6];

marker_size = 36; 
marker_alpha = 0.6; 

plot_suffix_no_outliers = '_TV_compare_scatter_jitter_sig'; 

% 显著性检验参数
alpha_level = 0.05; % 显著性水平
min_samples_for_test = 3; % 每组进行统计检验所需的最少样本数

for i = 1:length(region_names)
    current_region = region_names{i};
    fprintf('正在处理区域: %s\n', current_region);
    try
        neuron_data_var_name = ['TV_neuron_value_within_' current_region];
        edge_data_var_name = ['TV_value_within_' current_region];
        data_neuron_raw = eval(neuron_data_var_name);
        data_neuron_raw = data_neuron_raw(:); 
        data_edge_raw = eval(edge_data_var_name);
        data_edge_raw = data_edge_raw(:); 

        % --- 离群值移除 ---
        data_neuron_cleaned = data_neuron_raw;
        if ~isempty(data_neuron_raw)
            tf_neuron = isoutlier(data_neuron_raw, 'quartiles');
            data_neuron_cleaned = data_neuron_raw(~tf_neuron);
            fprintf('区域 %s, neuron 组: 原始数据点 %d 个, 移除离群值后 %d 个 (移除了 %d 个).\n', ...
                    current_region, length(data_neuron_raw), length(data_neuron_cleaned), sum(tf_neuron));
        end

        data_edge_cleaned = data_edge_raw;
        if ~isempty(data_edge_raw)
            tf_edge = isoutlier(data_edge_raw, 'quartiles');
            data_edge_cleaned = data_edge_raw(~tf_edge);
            fprintf('区域 %s, edge 组: 原始数据点 %d 个, 移除离群值后 %d 个 (移除了 %d 个).\n', ...
                    current_region, length(data_edge_raw), length(data_edge_cleaned), sum(tf_edge));
        end

        if isempty(data_neuron_cleaned) && isempty(data_edge_cleaned)
            fprintf('区域 %s 在移除离群值后，neuron 和 edge 数据均为空，跳过绘图。\n', current_region);
            continue;
        end

        figure;
        hold on;
        num_neuron_points = length(data_neuron_cleaned);
        if num_neuron_points > 0
            x_neuron = 1 + (rand(num_neuron_points, 1) - 0.5) * jitter_amount * 2;
            scatter(x_neuron, data_neuron_cleaned, marker_size, ...
                'filled', ...
                'MarkerFaceColor', neuron_color, ...
                'MarkerFaceAlpha', marker_alpha, ...
                'DisplayName', 'neuron'); 
        end

        num_edge_points = length(data_edge_cleaned);
        if num_edge_points > 0
            x_edge = 2 + (rand(num_edge_points, 1) - 0.5) * jitter_amount * 2;
            scatter(x_edge, data_edge_cleaned, marker_size, ...
                'filled', ...
                'MarkerFaceColor', edge_color, ...
                'MarkerFaceAlpha', marker_alpha, ...
                'DisplayName', 'edge'); 
        end

        perform_test_and_plot_sig = false; 
        if num_neuron_points >= min_samples_for_test && num_edge_points >= min_samples_for_test
            perform_test_and_plot_sig = true;
            % [h, p_value] = ttest2(data_neuron_cleaned, data_edge_cleaned, 'Vartype','unequal'); % 假设方差不相等
            [p_value, h_ranksum] = ranksum(data_neuron_cleaned, data_edge_cleaned); % h_ranksum 为 1 表示拒绝原假设 (即有显著差异)
            fprintf('区域 %s, Neuron vs Edge: p-value = %.4f (Wilcoxon rank-sum test)\n', current_region, p_value);
            % 确定显著性标记
            if p_value < 0.001
                sig_symbol = '***';
            elseif p_value < 0.01
                sig_symbol = '**';
            elseif p_value < alpha_level % alpha_level 通常是 0.05
                sig_symbol = '*';
            else
                sig_symbol = 'ns'; % non-significant
            end

            y_max_all_data = -Inf;
            if ~isempty(data_neuron_cleaned)
                y_max_all_data = max(y_max_all_data, max(data_neuron_cleaned));
            end
            if ~isempty(data_edge_cleaned)
                y_max_all_data = max(y_max_all_data, max(data_edge_cleaned));
            end

            if isfinite(y_max_all_data) 
                current_ylim = ylim; 
                y_range = current_ylim(2) - current_ylim(1);
                if y_range <= 0 || ~isfinite(y_range) 
                    y_range = 1; 
                end

                line_y_pos = y_max_all_data + 0.05 * y_range; 
                text_y_pos = line_y_pos + 0.03 * y_range;   

                % 绘制连接线
                plot([1, 2], [line_y_pos, line_y_pos], '-k', 'LineWidth', 1); 

                plot([1, 1], [line_y_pos - 0.01*y_range, line_y_pos], '-k', 'LineWidth', 1);
                plot([2, 2], [line_y_pos - 0.01*y_range, line_y_pos], '-k', 'LineWidth', 1);

                text(1.5, text_y_pos, sig_symbol, ...
                     'HorizontalAlignment', 'center', ...
                     'VerticalAlignment', 'bottom', ...
                     'FontSize', 12, 'FontWeight', 'bold');

                if text_y_pos > current_ylim(2)
                    ylim([current_ylim(1), text_y_pos + 0.05 * y_range]);
                end
            else
                 fprintf('区域 %s, 无法确定有效的 y_max_all_data 来绘制显著性标记。\n', current_region);
            end

        elseif num_neuron_points > 0 || num_edge_points > 0 
             fprintf('区域 %s, Neuron (%d) 或 Edge (%d) 组样本数不足 %d，跳过显著性检验。\n', ...
                     current_region, num_neuron_points, num_edge_points, min_samples_for_test);
        end

        hold off;

        title_str = [current_region ' region TV compare (Jitter, Outliers Removed)'];
        if perform_test_and_plot_sig
             title_str = [title_str, sprintf(', p=%.3f', p_value)]; 
        end
        title(title_str);
        ylabel('trial variability', 'FontSize', 12, 'FontName', 'Arial','FontWeight', 'bold');
        xlabel('data groups', 'FontSize', 12, 'FontName', 'Arial','FontWeight', 'bold');

        xlim([0.5, 2.5]); 
        xticks([1, 2]);   
        xticklabels({'SN', 'SNFC'}); 

        grid on; 
        box on;  

        file_basename = fullfile(figure_path, sprintf('%s%ssti%d', current_region, plot_suffix_no_outliers,show_std_sti_i));
        savefig([file_basename '.fig']);
        exportgraphics(gcf, [file_basename '.png'], 'Resolution', 300);
        fprintf('已保存 %s 的图形 (移除了离群值，并添加了显著性检验)。\n', current_region);

    catch ME
        fprintf('处理区域 %s 时发生错误: %s\n', current_region, ME.message);
        disp(ME.getReport());
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% jitter plot
% 定义脑区列表
region_names = {'VISp', 'VIS', 'VISa','VISal','VISam','VISl','VISli','VISpl','VISpor','VISrl','VISpm','SS', 'RSP', 'MO', 'AUD'};

jitter_amount = 0.15; % 控制抖动的幅度，可以根据数据点数量调整

neuron_color = [0.2 0.5 0.2];
edge_color = [1 0.85 0.6];

marker_size = 36; 
marker_alpha = 0.6; 

plot_suffix_no_outliers = '_TV_compare_scatter_jitter'; % 新的文件名后缀

for i = 1:length(region_names)
    current_region = region_names{i};
    fprintf('正在处理区域: %s\n', current_region);

    try
        neuron_data_var_name = ['TV_neuron_value_within_' current_region];
        edge_data_var_name = ['TV_value_within_' current_region];

        data_neuron_raw = eval(neuron_data_var_name);
        data_neuron_raw = data_neuron_raw(:); 

        data_edge_raw = eval(edge_data_var_name);
        data_edge_raw = data_edge_raw(:); % 确保是列向量

        % 对 neuron 数据移除离群值
        data_neuron_cleaned = data_neuron_raw; % 默认为原始数据
        if ~isempty(data_neuron_raw)
            tf_neuron = isoutlier(data_neuron_raw, 'quartiles');
            data_neuron_cleaned = data_neuron_raw(~tf_neuron);
            fprintf('区域 %s, neuron 组: 原始数据点 %d 个, 移除离群值后 %d 个 (移除了 %d 个).\n', ...
                    current_region, length(data_neuron_raw), length(data_neuron_cleaned), sum(tf_neuron));
        end

        % 对 edge 数据移除离群值
        data_edge_cleaned = data_edge_raw; 
        if ~isempty(data_edge_raw)
            tf_edge = isoutlier(data_edge_raw, 'quartiles');
            data_edge_cleaned = data_edge_raw(~tf_edge);
            fprintf('区域 %s, edge 组: 原始数据点 %d 个, 移除离群值后 %d 个 (移除了 %d 个).\n', ...
                    current_region, length(data_edge_raw), length(data_edge_cleaned), sum(tf_edge));
        end

        % 检查移除离群值后数据是否为空
        if isempty(data_neuron_cleaned) && isempty(data_edge_cleaned)
            fprintf('区域 %s 在移除离群值后，neuron 和 edge 数据均为空，跳过绘图。\n', current_region);
            continue;
        end

        figure;
        hold on;

        % 绘制 "neuron" 组数据 (使用清理后的数据)
        num_neuron_points = length(data_neuron_cleaned);
        if num_neuron_points > 0
            x_neuron = 1 + (rand(num_neuron_points, 1) - 0.5) * jitter_amount * 2;
            scatter(x_neuron, data_neuron_cleaned, marker_size, ...
                'filled', ...
                'MarkerFaceColor', neuron_color, ...
                'MarkerFaceAlpha', marker_alpha, ...
                'DisplayName', 'neuron ');
        end

        % 绘制 "edge" 组数据 (使用清理后的数据)
        num_edge_points = length(data_edge_cleaned);
        if num_edge_points > 0
            x_edge = 2 + (rand(num_edge_points, 1) - 0.5) * jitter_amount * 2;
            scatter(x_edge, data_edge_cleaned, marker_size, ...
                'filled', ...
                'MarkerFaceColor', edge_color, ...
                'MarkerFaceAlpha', marker_alpha, ...
                'DisplayName', 'edge ');
        end

        hold off;

        title([current_region ' region TV compare (Jitter Plot, Outliers Removed)']); % 更新标题
        ylabel('trial variability', 'FontSize', 12, 'FontName', 'Arial','FontWeight', 'bold');
        xlabel('data groups', 'FontSize', 12, 'FontName', 'Arial','FontWeight', 'bold');

        xlim([0.5, 2.5]); 
        xticks([1, 2]);   
        xticklabels({'SN', 'SNFC'}); 

        grid on; 
        box on;  

        file_basename = [figure_path sprintf('%s%s%sti%d', current_region, plot_suffix_no_outliers,show_std_sti_i)];
        savefig([file_basename '.fig']);
        exportgraphics(gcf, [file_basename '.png'], 'Resolution', 300);
        fprintf('已保存 %s 的图形 (移除了离群值)。\n', current_region);

    catch ME
        fprintf('处理区域 %s 时发生错误: %s\n', current_region, ME.message);
    end
end




% 定义脑区列表
% region_names = {'VISp', 'VIS', 'VISa','VISal','VISam','VISl','VISli','VISpl','VISpor','VISrl','VISpm','SS', 'RSP', 'MO', 'AUD'};
region_names = {'VISp', 'VISa','VISal','VISam','VISl','VISli','VISpl','VISpor','VISrl','VISpm','SS', 'RSP', 'MO', 'AUD'};

jitter_amount = 0.15;
neuron_color = [0.2 0.5 0.2];
edge_color = [1 0.85 0.6];
marker_size = 36;
marker_alpha = 0.6;
alpha_level = 0.05;
min_samples_for_test = 3;

all_neuron_data = {};
all_edge_data = {};
all_regions_used = {};
sig_ps = [];

% 先将符合要求的区域数据选出来
for i = 1:length(region_names)
    current_region = region_names{i};
    neuron_data_var_name = ['TV_neuron_value_within_' current_region];
    edge_data_var_name = ['TV_value_within_' current_region];

    % 获取数据并离群值处理
    data_neuron = eval(neuron_data_var_name);
    data_edge = eval(edge_data_var_name);
    data_neuron = data_neuron(:);
    data_edge   = data_edge(:);

    if ~isempty(data_neuron)
        data_neuron = data_neuron(~isoutlier(data_neuron, 'quartiles'));
    end
    if ~isempty(data_edge)
        data_edge = data_edge(~isoutlier(data_edge, 'quartiles'));
    end

    % 样本数量要求
    if (length(data_neuron) >= min_samples_for_test) && (length(data_edge) >= min_samples_for_test)
        % 显著性检验
        [p, ~] = ranksum(data_neuron, data_edge);
        all_neuron_data{end+1} = data_neuron;
        all_edge_data{end+1} = data_edge;
        all_regions_used{end+1} = current_region;
        sig_ps(end+1) = p; %#ok<SAGROW>
    end
end


figure; hold on;
xticklabels_cell = {};
bar_locs = [];
for i = 1:length(all_regions_used)
    neuron_data = all_neuron_data{i};
    edge_data   = all_edge_data{i};
    region_name = all_regions_used{i};

    % --- 为了防止数据中出现0或负数导致对数坐标轴错误，过滤掉非正值 ---
    neuron_data(neuron_data <= 0) = NaN;
    edge_data(edge_data <= 0) = NaN;

    x_neuron = i*2-1 + (rand(length(neuron_data),1)-0.5)*jitter_amount*2;
    x_edge   = i*2   + (rand(length(edge_data),1)-0.5)*jitter_amount*2;
    scatter(x_neuron, neuron_data, marker_size, 'filled', 'MarkerFaceColor', neuron_color, 'MarkerFaceAlpha', marker_alpha);
    scatter(x_edge  , edge_data  , marker_size, 'filled', 'MarkerFaceColor', edge_color , 'MarkerFaceAlpha', marker_alpha);

    bar_locs = [bar_locs, i*2-1, i*2];
    xticklabels_cell = [xticklabels_cell, [region_name '-SN'], [region_name '-SNFC']];
end

set(gca, 'YScale', 'log');

% 显著性线
for i = 1:length(all_regions_used)
    neuron_data = all_neuron_data{i};
    edge_data = all_edge_data{i};
    p = sig_ps(i);

    if p < 0.001
        sig_symbol = '***';
    elseif p < 0.01
        sig_symbol = '**';
    elseif p < alpha_level
        sig_symbol = '*';
    else
        sig_symbol = 'ns';
    end

    ymax = max([neuron_data; edge_data]);

    % 如果当前组的最大值为NaN（例如所有数据都被过滤），则跳过该组的显著性标记
    if isnan(ymax)
        continue;
    end   

    line_factor = 1.2;  
    text_factor = 1.1;  
    tick_factor = 0.98; 

    line_y = ymax * line_factor;
    text_y = line_y * text_factor;

    x1 = i*2-1; x2 = i*2;
    plot([x1, x2], [line_y, line_y], '-k', 'LineWidth', 1);
    plot([x1, x1], [line_y * tick_factor, line_y], '-k', 'LineWidth', 1); 
    plot([x2, x2], [line_y * tick_factor, line_y], '-k', 'LineWidth', 1); 
    text((x1+x2)/2, text_y, sig_symbol, 'HorizontalAlignment','center','FontSize',12,'FontWeight','bold');

    cur_ylim = ylim;
    if text_y > cur_ylim(2)
        ylim([cur_ylim(1), text_y * 1.2]); 
    end
end

% ======================================================

set(gca, 'XTick', bar_locs, 'XTickLabel', xticklabels_cell, 'XTickLabelRotation', 45);
ylabel('trial variability (log scale)','FontSize',12,'FontWeight','bold');
title('Significant Regions TV Compare (Jitter, Outliers Removed)');
grid on; box on; hold off;
set(gcf,'Position',[100 100 1400 600]);

savefig(fullfile(figure_path,sprintf('Significant_Regions_TV_Compare_jitter_sig_sti%d_log.fig',show_std_sti_i)));
exportgraphics(gcf, fullfile(figure_path,sprintf('Significant_Regions_TV_Compare_jitter_sig_sti%d_log.png',show_std_sti_i)),'Resolution',300);


close all

end




















%% ----------------------------------------------------------------------------------------------------------------------------------------------- show corr matrix for similar select results
% multi trace average imagesc
[C_mat_cell,before,during,after,sti_num] = get_trials_matrix_multi_trace_by_sti((C_mat),stimuli,param);

param_for_corr = param;
param_for_corr.show_sti_before = 6;
[select_neurons_edge_cell,before_for_corr,during_for_corr,after_for_corr,sti_num] = get_trials_matrix_multi_trace_by_sti(select_neurons_C_sort_by_region_corr_matrix,stimuli,param_for_corr);

% >>>>>>>>>>>>>>>>>>> >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> sort average trial (neuron)
% >>>>>>>>>>>>>>>>>>> sort average trial all neurons
close all
average_trial_sti_select_neurons_cell = {};
for sti_i = 1:size(sort_all_neuron_stdResult_id,2)
    curr_select_id_buff = sort_all_neuron_stdResult_id{1,sti_i};

    curr_select_id = curr_select_id_buff(1:select_neuron_num,:);
    figure;
    clim_val = [-0.1 0.2]; 
    for data_sessioni = 1:data_session_n
        % curr_sti_select_neurons_cell = C_matNorm_cell{1,sti_i};
        curr_sti_select_neurons_cell = C_mat_cell{1,sti_i};
        average_trial_sti_select_neurons_cell_buff = zeros(size(curr_sti_select_neurons_cell{1,1}));

        for trial_i = each_session_trials_num*(data_sessioni-1)+1:each_session_trials_num*data_sessioni
            curr_trial_sti_select_neurons_cell = curr_sti_select_neurons_cell{1,trial_i};
            average_trial_sti_select_neurons_cell_buff = average_trial_sti_select_neurons_cell_buff+curr_trial_sti_select_neurons_cell;
        end

        average_trial_sti_select_neurons_cell0 = average_trial_sti_select_neurons_cell_buff/each_session_trials_num;
        average_trial_sti_select_neurons_cell{sti_i,data_sessioni} = average_trial_sti_select_neurons_cell0;
        average_trial_sti_select_neurons_select_cell0 = average_trial_sti_select_neurons_cell0(curr_select_id,:);
        average_trial_sti_select_neurons_select_cell{sti_i,data_sessioni} = average_trial_sti_select_neurons_select_cell0;


        average_trial_sti_select_neurons_select_cell_first = average_trial_sti_select_neurons_select_cell{sti_i,1}(:,before+1:end);
        [average_trial_sti_select_neurons_select_cell_sort,select_average_sortedTimes, select_average_sortedValues,select_average_sortIdx] = ...
            sort_trace_by_peakTime(average_trial_sti_select_neurons_select_cell_first');    
        average_trial_sti_select_neurons_select_cell0_sort = average_trial_sti_select_neurons_select_cell0(select_average_sortIdx,:);
        subplot(1,data_session_n,data_sessioni);

        imagesc(zscore(average_trial_sti_select_neurons_select_cell0_sort, 0, 2))
        colormap(red_blue_color);
        % clim(clim_val);
        y_axis  = 0.5 : 0.5 : size(average_trial_sti_select_neurons_select_cell0_sort, 1);
        hold on, plot(before * ones(1, length(y_axis)), y_axis, '--r', 'linewidth', 2)
        hold on, plot((before+during) * ones(1, length(y_axis)), y_axis, '--r', 'linewidth', 2)

        xticks = get(gca, 'XTick'); 
        xticklabels = xticks / param.fs - param.show_sti_before; 
        set(gca, 'XTick', xticks, 'XTickLabel', xticklabels); 


        title(sprintf('sti %d',sti_i), 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold'); 
        xlabel('Time (s)', 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');    
        ylabel('Neuron #', 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');   
        set(gca, 'FontSize', 8, 'FontName', 'Arial');            

    end
    c = colorbar('Position', [0.92 0.1 0.02 0.82]); 
    ylabel(c, 'neuron activity intensity value','FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold'); 

    top_num = size(average_trial_sti_select_neurons_select_cell0_sort,1);
    savefig([figure_path sprintf('top_%d_stable_neuron_sort_sti_%d_average_trial',top_num,sti_i)]);
    exportgraphics(gcf,[figure_path sprintf('top_%d_stable_neuron_sort_sti_%d_average_trial',top_num,sti_i),'.png'],'Resolution',300)
end


% (plot in one fig)
close all
average_trial_sti_select_neurons_cell = {};
figure;
clim_val = [-0.1 0.2]; 
for sti_i = 1:size(sort_all_neuron_stdResult_id,2)
    % curr_select_id_buff = sort_all_neuron_stdResult_id{1,sti_i};
    curr_select_id_buff = sort_all_neuron_stdResult_id{1,1};
    curr_select_id = curr_select_id_buff(1:select_neuron_num,:);

    % curr_sti_select_neurons_cell = C_matNorm_cell{1,sti_i};
    curr_sti_select_neurons_cell = C_mat_cell{1,sti_i};
    average_trial_sti_select_neurons_cell_buff = zeros(size(curr_sti_select_neurons_cell{1,1}));

    for trial_i = each_session_trials_num*(1-1)+1:each_session_trials_num*1
        curr_trial_sti_select_neurons_cell = curr_sti_select_neurons_cell{1,trial_i};
        average_trial_sti_select_neurons_cell_buff = average_trial_sti_select_neurons_cell_buff+curr_trial_sti_select_neurons_cell;
    end

    average_trial_sti_select_neurons_cell0 = average_trial_sti_select_neurons_cell_buff/each_session_trials_num;
    average_trial_sti_select_neurons_cell{sti_i,1} = average_trial_sti_select_neurons_cell0;
    average_trial_sti_select_neurons_select_cell0 = average_trial_sti_select_neurons_cell0(curr_select_id,:);
    average_trial_sti_select_neurons_select_cell{sti_i,data_sessioni} = average_trial_sti_select_neurons_select_cell0;

    average_trial_sti_select_neurons_select_cell_first = average_trial_sti_select_neurons_select_cell{sti_i,1}(:,before+1:end);
    [average_trial_sti_select_neurons_select_cell_sort,select_average_sortedTimes, select_average_sortedValues,select_average_sortIdx] = ...
        sort_trace_by_peakTime(average_trial_sti_select_neurons_select_cell_first');    
    average_trial_sti_select_neurons_select_cell0_sort = average_trial_sti_select_neurons_select_cell0(select_average_sortIdx,:);
    subplot(1,size(sort_all_neuron_stdResult_id,2),sti_i);

    % imagesc(average_trial_sti_select_neurons_select_cell0_sort); hold on;
    imagesc(zscore(average_trial_sti_select_neurons_select_cell0_sort, 0, 2));hold on;
    colormap(red_blue_color);
    % clim(clim_val);
    y_axis  = 0.5 : 0.5 : size(average_trial_sti_select_neurons_select_cell0_sort, 1);
    hold on, plot(before * ones(1, length(y_axis)), y_axis, '--r', 'linewidth', 2)
    hold on, plot((before+during) * ones(1, length(y_axis)), y_axis, '--r', 'linewidth', 2)

    xticks = get(gca, 'XTick'); 
    xticklabels = xticks / param.fs - param.show_sti_before; 
    set(gca, 'XTick', xticks, 'XTickLabel', xticklabels); 

    title(sprintf('sti %d',sti_i), 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold'); 
    xlabel('Time (s)', 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');    
    ylabel('Neuron #', 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');    
    set(gca, 'FontSize', 8, 'FontName', 'Arial');            

end
c = colorbar('Position', [0.92 0.1 0.02 0.82]); 
ylabel(c, 'neuron activity intensity value','FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold'); 

top_num = size(average_trial_sti_select_neurons_select_cell0_sort,1);
savefig([figure_path sprintf('top_%d_stable_neuron_sortbysti1_all_sti_average_trial',top_num)]);
exportgraphics(gcf,[figure_path sprintf('top_%d_stable_neuron_sortbysti1_all_sti_average_trial',top_num),'.png'],'Resolution',300)

close all
% >>>>>>>>>>>>>>>>>>> sort each trial select neurons 
for sti_i = 1:size(sort_all_neuron_stdResult_id,2)
    curr_select_id_buff = sort_all_neuron_stdResult_id{1,sti_i};
    curr_select_id = curr_select_id_buff(1:select_neuron_num,:);
    % curr_sti_select_neurons_cell = C_matNorm_cell{1,sti_i};
    curr_sti_select_neurons_cell = C_mat_cell{1,sti_i};
    % sort trace by first trial (annotating this means use average sort)
    first_trial_sti_select_neurons_cell = curr_sti_select_neurons_cell{1,1}(curr_select_id,before+1:end);
    [first_trial_sti_select_neurons_cell_sort,first_sortedTimes, first_sortedValues,first_sortIdx] = sort_trace_by_peakTime(first_trial_sti_select_neurons_cell');

    figure;
    clim_val = [-0.1 0.2];
    step = 4;
    n = 0;
    % for trial_i = 1:step:size(curr_sti_select_neurons_cell,2)
    for trial_i = 1:5
        n = n+1;
        curr_trial_sti_select_neurons_cell = curr_sti_select_neurons_cell{1,trial_i}(curr_select_id,:);
        curr_trial_sti_select_neurons_cell_sort = curr_trial_sti_select_neurons_cell(first_sortIdx,:);
        % curr_trial_sti_select_neurons_cell_sort = curr_trial_sti_select_neurons_cell(average_sortIdx,:);
        % subplot(1,size(curr_sti_select_neurons_cell,2)/step,n);
        subplot(1,5,n);
        imagesc(zscore(curr_trial_sti_select_neurons_cell_sort, 0, 2));hold on;
        colormap(red_blue_color);
        % clim(clim_val);
        y_axis  = 0.5 : 0.5 : size(curr_trial_sti_select_neurons_cell_sort, 1);
        hold on, plot(before * ones(1, length(y_axis)), y_axis, '--r', 'linewidth', 2)
        hold on, plot((before+during) * ones(1, length(y_axis)), y_axis, '--r', 'linewidth', 2)

        xticks = get(gca, 'XTick'); 
        xticklabels = xticks / param.fs - param.show_sti_before; 
        set(gca, 'XTick', xticks, 'XTickLabel', xticklabels); 

        title(sprintf('trial %d',trial_i), 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');
        xlabel('Time (s)', 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');  
        ylabel('Neuron #', 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');   
        set(gca, 'FontSize', 8, 'FontName', 'Arial');           
    end
    c = colorbar('Position', [0.92 0.1 0.02 0.82]); 
    ylabel(c, 'neuron activity intensity value','FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold'); 
    top_num = size(curr_trial_sti_select_neurons_cell_sort,1);
    savefig([figure_path sprintf('multi_session_select_%d_neuron_sort_sti_%d_1to5single_trial',top_num,sti_i)]);
    exportgraphics(gcf,[figure_path sprintf('multi_session_select_%d_neuron_sort_sti_%d_1to5_single_trial',top_num,sti_i),'.png'],'Resolution',300)
end


% select edge (select by ttest result) 
for sti_i = 1:size(sort_select_neuron_stdResult_id,2)
    curr_select_id_buff = sort_select_neuron_stdResult_id{1,sti_i};
    curr_select_id = curr_select_id_buff(1:select_neuron_num,:);
    figure;
    clim_val = [-0.1 0.3];
    for data_sessioni = 1:data_session_n
        % curr_sti_select_neurons_cell = C_matNorm_cell{1,sti_i};
        curr_sti_select_neurons_cell = C_mat_cell{1,sti_i};
        average_trial_sti_select_neurons_cell_buff = zeros(size(curr_sti_select_neurons_cell{1,1}));

        for trial_i = each_session_trials_num*(data_sessioni-1)+1:each_session_trials_num*data_sessioni
            curr_trial_sti_select_neurons_cell = curr_sti_select_neurons_cell{1,trial_i};
            average_trial_sti_select_neurons_cell_buff = average_trial_sti_select_neurons_cell_buff+curr_trial_sti_select_neurons_cell;
        end

        average_trial_sti_select_neurons_cell0 = average_trial_sti_select_neurons_cell_buff/each_session_trials_num;
        average_trial_sti_select_neurons_cell{sti_i,data_sessioni} = average_trial_sti_select_neurons_cell0;
        average_trial_sti_select_neurons_select_cell0 = average_trial_sti_select_neurons_cell0(select_neurons_id(curr_select_id),:);
        average_trial_sti_select_neurons_select_cell{sti_i,data_sessioni} = average_trial_sti_select_neurons_select_cell0;

        average_trial_sti_select_neurons_select_cell_first = average_trial_sti_select_neurons_select_cell{sti_i,1}(:,before+1:end);
        [average_trial_sti_select_neurons_select_cell_sort,select_average_sortedTimes, select_average_sortedValues,select_average_sortIdx] = ...
            sort_trace_by_peakTime(average_trial_sti_select_neurons_select_cell_first');    
        average_trial_sti_select_neurons_select_cell0_sort = average_trial_sti_select_neurons_select_cell0(select_average_sortIdx,:);
        subplot(1,data_session_n,data_sessioni);
        imagesc(zscore(average_trial_sti_select_neurons_select_cell0_sort, 0, 2));hold on;
        colormap(red_blue_color);
        % clim(clim_val);
        y_axis  = 0.5 : 0.5 : size(average_trial_sti_select_neurons_select_cell0_sort, 1);
        hold on, plot(before * ones(1, length(y_axis)), y_axis, '--r', 'linewidth', 2)
        hold on, plot((before+during) * ones(1, length(y_axis)), y_axis, '--r', 'linewidth', 2)

        xticks = get(gca, 'XTick'); 
        xticklabels = xticks / param.fs - param.show_sti_before; 
        set(gca, 'XTick', xticks, 'XTickLabel', xticklabels); 

        title(sprintf('sti %d',sti_i), 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold'); 
        xlabel('Time (s)', 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');    
        ylabel('Neuron #', 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');   
        set(gca, 'FontSize', 8, 'FontName', 'Arial');         
    end
    c = colorbar('Position', [0.92 0.1 0.02 0.82]); 
    ylabel(c, 'neuron activity intensity value','FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold'); 
    top_num = size(average_trial_sti_select_neurons_select_cell0_sort,1);
    savefig([figure_path sprintf('multi_session_top_%d_ttest_stable_neuron_sort_sti_%d_average_trial',top_num,sti_i)]);
    exportgraphics(gcf,[figure_path sprintf('multi_session_top_%d_ttest_stable_neuron_sort_sti_%d_average_trial',top_num,sti_i),'.png'],'Resolution',300)
end

close all
figure;
clim_val = [-0.1 0.3];
% select edge (select by ttest result)
for sti_i = 1:size(sort_select_neuron_stdResult_id,2)
    % curr_select_id_buff = sort_select_neuron_stdResult_id{1,sti_i};
    curr_select_id_buff = sort_select_neuron_stdResult_id{1,1};
    curr_select_id = curr_select_id_buff(1:select_neuron_num,:);
    data_sessioni = 1;
        % curr_sti_select_neurons_cell = C_matNorm_cell{1,sti_i};
        curr_sti_select_neurons_cell = C_mat_cell{1,sti_i};
        average_trial_sti_select_neurons_cell_buff = zeros(size(curr_sti_select_neurons_cell{1,1}));

        for trial_i = each_session_trials_num*(data_sessioni-1)+1:each_session_trials_num*data_sessioni
            curr_trial_sti_select_neurons_cell = curr_sti_select_neurons_cell{1,trial_i};
            average_trial_sti_select_neurons_cell_buff = average_trial_sti_select_neurons_cell_buff+curr_trial_sti_select_neurons_cell;
        end

        average_trial_sti_select_neurons_cell0 = average_trial_sti_select_neurons_cell_buff/each_session_trials_num;
        average_trial_sti_select_neurons_cell{sti_i,data_sessioni} = average_trial_sti_select_neurons_cell0;
        average_trial_sti_select_neurons_select_cell0 = average_trial_sti_select_neurons_cell0(select_neurons_id(curr_select_id),:);
        average_trial_sti_select_neurons_select_cell{sti_i,data_sessioni} = average_trial_sti_select_neurons_select_cell0;

        average_trial_sti_select_neurons_select_cell_first = average_trial_sti_select_neurons_select_cell{sti_i,1}(:,before+1:end);
        [average_trial_sti_select_neurons_select_cell_sort,select_average_sortedTimes, select_average_sortedValues,select_average_sortIdx] = ...
            sort_trace_by_peakTime(average_trial_sti_select_neurons_select_cell_first');    
        average_trial_sti_select_neurons_select_cell0_sort = average_trial_sti_select_neurons_select_cell0(select_average_sortIdx,:);
        subplot(1,size(sort_select_neuron_stdResult_id,2),sti_i);
        % imagesc(average_trial_sti_select_neurons_select_cell0_sort)
        imagesc(zscore(average_trial_sti_select_neurons_select_cell0_sort, 0, 2));hold on;
        colormap(red_blue_color);
        % clim(clim_val);
        y_axis  = 0.5 : 0.5 : size(average_trial_sti_select_neurons_select_cell0_sort, 1);
        hold on, plot(before * ones(1, length(y_axis)), y_axis, '--r', 'linewidth', 2)
        hold on, plot((before+during) * ones(1, length(y_axis)), y_axis, '--r', 'linewidth', 2)

        xticks = get(gca, 'XTick'); 
        xticklabels = xticks / param.fs - param.show_sti_before; 
        set(gca, 'XTick', xticks, 'XTickLabel', xticklabels); 

        title(sprintf('sti %d',sti_i), 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold'); 
        xlabel('Time (s)', 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');    
        ylabel('Neuron #', 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');    
        set(gca, 'FontSize', 8, 'FontName', 'Arial');          

end
c = colorbar('Position', [0.92 0.1 0.02 0.82]); 
ylabel(c, 'neuron activity intensity value','FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold'); 
top_num = size(average_trial_sti_select_neurons_select_cell0_sort,1);
savefig([figure_path sprintf('top_%d_ttest_stable_neuron_sort_all_sti_average_trial',top_num)]);
exportgraphics(gcf,[figure_path sprintf('top_%d_ttest_stable_neuron_sort_all_sti_average_trial',top_num),'.png'],'Resolution',300)

close all
% >>>>>>>>>>>>>>>>>>> sort each trial select neurons
for sti_i = 1:size(sort_all_neuron_stdResult_id,2)
    curr_select_id_buff = sort_select_neuron_stdResult_id{1,sti_i};
    curr_select_id = curr_select_id_buff(1:select_neuron_num,:);
    % curr_select_id = curr_select_id_buff;
    % curr_sti_select_neurons_cell = C_matNorm_cell{1,sti_i};
    curr_sti_select_neurons_cell = C_mat_cell{1,sti_i};
    % sort trace by first trial (annotating this means use average sort)
    first_trial_sti_select_neurons_cell = curr_sti_select_neurons_cell{1,1}(curr_select_id,before+1:end);
    [first_trial_sti_select_neurons_cell_sort,first_sortedTimes, first_sortedValues,first_sortIdx] = sort_trace_by_peakTime(first_trial_sti_select_neurons_cell');
    figure;
    clim_val = [-0.05 0.3]; 
    step = 4;
    n = 0;
    % for trial_i = 1:step:size(curr_sti_select_neurons_cell,2)
    for trial_i = 1:5
        n = n+1;
        curr_trial_sti_select_neurons_cell = curr_sti_select_neurons_cell{1,trial_i}(select_neurons_id(curr_select_id),:);
        curr_trial_sti_select_neurons_cell_sort = curr_trial_sti_select_neurons_cell(first_sortIdx,:);
        % curr_trial_sti_select_neurons_cell_sort = curr_trial_sti_select_neurons_cell(average_sortIdx,:);
        subplot(1,5,n);
        imagesc(zscore(curr_trial_sti_select_neurons_cell_sort, 0, 2));hold on;
        colormap(red_blue_color);
        % clim(clim_val);
        y_axis  = 0.5 : 0.5 : size(curr_trial_sti_select_neurons_cell_sort, 1);
        hold on, plot(before * ones(1, length(y_axis)), y_axis, '--r', 'linewidth', 2)
        hold on, plot((before+during) * ones(1, length(y_axis)), y_axis, '--r', 'linewidth', 2)

        xticks = get(gca, 'XTick'); 
        xticklabels = xticks / param.fs - param.show_sti_before; 
        set(gca, 'XTick', xticks, 'XTickLabel', xticklabels); 

        title(sprintf('trial %d',trial_i), 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold'); 
        xlabel('Time (s)', 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');    
        ylabel('Neuron #', 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');    
        set(gca, 'FontSize', 8, 'FontName', 'Arial');           

    end
    c = colorbar('Position', [0.92 0.1 0.02 0.82]); 
    ylabel(c, 'neuron activity intensity value','FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold'); 
    top_num = size(curr_trial_sti_select_neurons_cell_sort,1);
    savefig([figure_path sprintf('multi_session_select_%d_ttest_neuron_sort_sti_%d_1to5_single_trial',top_num,sti_i)]);
    exportgraphics(gcf,[figure_path sprintf('multi_session_select_%d_ttest_neuron_sort_sti_%d_1to5_single_trial',top_num,sti_i),'.png'],'Resolution',300)
end

close all
% >>>>>>>>>>>>>>>>>>> sort each trial all neurons
for sti_i = 1:size(sort_all_neuron_stdResult_id,2)
    % curr_sti_select_neurons_cell = C_matNorm_cell{1,sti_i};
    curr_sti_select_neurons_cell = C_mat_cell{1,sti_i};
    % sort trace by first trial (annotating this means use average sort)
    first_trial_sti_select_neurons_cell = curr_sti_select_neurons_cell{1,1}(:,before+1:end);
    [first_trial_sti_select_neurons_cell_sort,first_sortedTimes, first_sortedValues,first_sortIdx] = sort_trace_by_peakTime(first_trial_sti_select_neurons_cell');
    % figure;imagesc(first_trial_sti_select_neurons_cell_sort')
    % close all
    figure;
    % clim_val = [-100 200]; 
    step = 4;
    n = 0;
    % for trial_i = 1:step:size(curr_sti_select_neurons_cell,2)
    for trial_i = 1:5
        n = n+1;
        curr_trial_sti_select_neurons_cell = curr_sti_select_neurons_cell{1,trial_i};
        curr_trial_sti_select_neurons_cell_sort = curr_trial_sti_select_neurons_cell(first_sortIdx,:);
        % curr_trial_sti_select_neurons_cell_sort = curr_trial_sti_select_neurons_cell(average_sortIdx,:);
        % subplot(1,size(curr_sti_select_neurons_cell,2)/step,n);
        subplot(1,5,n);
        % imagesc(curr_trial_sti_select_neurons_cell_sort)
        imagesc(zscore(curr_trial_sti_select_neurons_cell_sort, 0, 2));hold on;
        colormap(red_blue_color);
        % clim(clim_val);
        y_axis  = 0.5 : 0.5 : size(curr_trial_sti_select_neurons_cell_sort, 1);
        hold on, plot(before * ones(1, length(y_axis)), y_axis, '--r', 'linewidth', 2)
        hold on, plot((before+during) * ones(1, length(y_axis)), y_axis, '--r', 'linewidth', 2)

        xticks = get(gca, 'XTick'); 
        xticklabels = xticks / param.fs - param.show_sti_before; 
        set(gca, 'XTick', xticks, 'XTickLabel', xticklabels); 

        title(sprintf('trial %d',trial_i), 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');
        xlabel('Time (s)', 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');   
        ylabel('Neuron #', 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');    
        set(gca, 'FontSize', 8, 'FontName', 'Arial');           
    end
    c = colorbar('Position', [0.92 0.1 0.02 0.82]); 
    ylabel(c, 'neuron activity intensity value','FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold'); 
    top_num = size(curr_trial_sti_select_neurons_cell_sort,1);
    savefig([figure_path sprintf('multi_session_all_%d_neuron_sort_sti_%d_1to5_single_trial',top_num,sti_i)]);
    exportgraphics(gcf,[figure_path sprintf('multi_session_all_%d_neuron_sort_sti_%d_1to5_single_trial',top_num,sti_i),'.png'],'Resolution',300)
end


% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> show edge

% >>>>>>>>>>>>>>>>>>> sort average trial all edge
close all
average_trial_sti_select_neurons_edge_cell = {};
for sti_i = 1:size(select_neurons_edge_cell,2)
    curr_edge_select_id_buff = sort_edge_stdResult_id{1,sti_i};
    curr_edge_select_id = curr_edge_select_id_buff(1:select_neuron_num,:);
    figure;
    clim_val = [-0.6 0.8]; 
    for data_sessioni = 1:data_session_n
        curr_sti_select_neurons_edge_cell = select_neurons_edge_cell{1,sti_i};
        average_trial_sti_select_neurons_edge_cell_buff = zeros(size(curr_sti_select_neurons_edge_cell{1,1}));

        for trial_i = each_session_trials_num*(data_sessioni-1)+1:each_session_trials_num*data_sessioni
            curr_trial_sti_select_neurons_edge_cell = curr_sti_select_neurons_edge_cell{1,trial_i};
            average_trial_sti_select_neurons_edge_cell_buff = average_trial_sti_select_neurons_edge_cell_buff+curr_trial_sti_select_neurons_edge_cell;
        end

        average_trial_sti_select_neurons_edge_cell0 = average_trial_sti_select_neurons_edge_cell_buff/each_session_trials_num;
        average_trial_sti_select_neurons_edge_cell{sti_i,data_sessioni} = average_trial_sti_select_neurons_edge_cell0;
        average_trial_sti_select_neurons_select_edge_cell0 = average_trial_sti_select_neurons_edge_cell0(curr_edge_select_id,:);
        average_trial_sti_select_neurons_select_edge_cell{sti_i,data_sessioni} = average_trial_sti_select_neurons_select_edge_cell0;

        average_trial_sti_select_neurons_select_edge_cell_first = average_trial_sti_select_neurons_select_edge_cell{sti_i,1};
        [average_trial_sti_select_neurons_select_edge_cell_sort,select_edge_average_sortedTimes, select_edge_average_sortedValues,select_edge_average_sortIdx] = ...
            sort_trace_by_peakTime(average_trial_sti_select_neurons_select_edge_cell_first');    
        average_trial_sti_select_neurons_select_edge_cell0_sort = average_trial_sti_select_neurons_select_edge_cell0(select_edge_average_sortIdx,:);
        subplot(1,data_session_n,data_sessioni);
        imagesc(zscore(average_trial_sti_select_neurons_select_edge_cell0_sort, 0, 2));hold on;
        colormap(red_blue_color);

        y_axis  = 0.5 : 0.5 : size(average_trial_sti_select_neurons_select_edge_cell0_sort, 1);
        hold on, plot(12 * ones(12, length(y_axis)), y_axis, '--r', 'linewidth', 2)
        hold on, plot(before_for_corr * ones(1, length(y_axis)), y_axis, '--w', 'linewidth', 2)
        hold on, plot((before_for_corr+during) * ones(1, length(y_axis)), y_axis, '--r', 'linewidth', 2)

        xticks = get(gca, 'XTick'); 
        xticklabels = xticks / param.fs - param.show_sti_before; 
        set(gca, 'XTick', xticks, 'XTickLabel', xticklabels); 

        title(sprintf('sti %d',sti_i), 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold'); 
        xlabel('Time (s)', 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');   
        ylabel('Edge #', 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');   
        set(gca, 'FontSize', 8, 'FontName', 'Arial');            
    end

    c = colorbar('Position', [0.92 0.1 0.02 0.82]); 
    ylabel(c, 'neuron activity intensity value','FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold'); 
    top_num = size(average_trial_sti_select_neurons_select_edge_cell0_sort,1);
    savefig([figure_path sprintf('multi_session_top_%d_stable_edge_sort_average_trial_sti_%d',top_num,sti_i)]);
    exportgraphics(gcf,[figure_path sprintf('multi_session_top_%d_stable_edge_sort_average_trial_sti_%d',top_num,sti_i),'.png'],'Resolution',300)
end

% >>>>>>>>>>>>>>>>>>> sort average trial all edge (plot in one fig)
close all
figure;
clim_val = [-0.6 0.8]; 
average_trial_sti_select_neurons_edge_cell = {};
for sti_i = 1:size(select_neurons_edge_cell,2)
    % curr_edge_select_id_buff = sort_edge_stdResult_id{1,sti_i};
    curr_edge_select_id_buff = sort_edge_stdResult_id{1,1};
    curr_edge_select_id = curr_edge_select_id_buff(1:select_neuron_num,:);

    % for data_sessioni = 1:data_session_n
    data_sessioni = 1;
    curr_sti_select_neurons_edge_cell = select_neurons_edge_cell{1,sti_i};
    average_trial_sti_select_neurons_edge_cell_buff = zeros(size(curr_sti_select_neurons_edge_cell{1,1}));

    for trial_i = each_session_trials_num*(data_sessioni-1)+1:each_session_trials_num*data_sessioni
        curr_trial_sti_select_neurons_edge_cell = curr_sti_select_neurons_edge_cell{1,trial_i};
        average_trial_sti_select_neurons_edge_cell_buff = average_trial_sti_select_neurons_edge_cell_buff+curr_trial_sti_select_neurons_edge_cell;
    end

    average_trial_sti_select_neurons_edge_cell0 = average_trial_sti_select_neurons_edge_cell_buff/each_session_trials_num;
    average_trial_sti_select_neurons_edge_cell{sti_i,data_sessioni} = average_trial_sti_select_neurons_edge_cell0;
    average_trial_sti_select_neurons_select_edge_cell0 = average_trial_sti_select_neurons_edge_cell0(curr_edge_select_id,:);
    average_trial_sti_select_neurons_select_edge_cell{sti_i,data_sessioni} = average_trial_sti_select_neurons_select_edge_cell0;

    % average_trial_sti_select_neurons_select_edge_cell_first = average_trial_sti_select_neurons_select_edge_cell{sti_i,1}(:,before+1:end);
    average_trial_sti_select_neurons_select_edge_cell_first = average_trial_sti_select_neurons_select_edge_cell{sti_i,1};
    [average_trial_sti_select_neurons_select_edge_cell_sort,select_edge_average_sortedTimes, select_edge_average_sortedValues,select_edge_average_sortIdx] = ...
        sort_trace_by_peakTime(average_trial_sti_select_neurons_select_edge_cell_first');    
    average_trial_sti_select_neurons_select_edge_cell0_sort = average_trial_sti_select_neurons_select_edge_cell0(select_edge_average_sortIdx,:);
    subplot(1,size(select_neurons_edge_cell,2),sti_i);
    imagesc(zscore(average_trial_sti_select_neurons_select_edge_cell0_sort, 0, 2));hold on;
    colormap(red_blue_color);
    y_axis  = 0.5 : 0.5 : size(average_trial_sti_select_neurons_select_edge_cell0_sort, 1);
    hold on, plot(12 * ones(12, length(y_axis)), y_axis, '--r', 'linewidth', 2)
    hold on, plot(before_for_corr * ones(1, length(y_axis)), y_axis, '--w', 'linewidth', 2)
    hold on, plot((before_for_corr+during) * ones(1, length(y_axis)), y_axis, '--r', 'linewidth', 2)

    xticks = get(gca, 'XTick'); 
    xticklabels = xticks / param.fs - param.show_sti_before; 
    set(gca, 'XTick', xticks, 'XTickLabel', xticklabels); 

    title(sprintf('sti %d',sti_i), 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold'); 
    xlabel('Time (s)', 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');    
    ylabel('Edge #', 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');   
    set(gca, 'FontSize', 8, 'FontName', 'Arial');            

end
c = colorbar('Position', [0.92 0.1 0.02 0.82]); 
ylabel(c, 'neuron activity intensity value','FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold'); 
top_num = size(average_trial_sti_select_neurons_select_edge_cell0_sort,1);
savefig([figure_path sprintf('multi_session_top_%d_stable_edge_sortbysti1_average_trial_all_sti',top_num)]);
exportgraphics(gcf,[figure_path sprintf('multi_session_top_%d_stable_edge_sortbysti1_average_trial_all_sti',top_num),'.png'],'Resolution',300)


close all
% >>>>>>>>>>>>>>>>>>> sort each single trial top edge
for sti_i = 1:size(select_neurons_edge_cell,2)

    curr_edge_select_id_buff = sort_edge_stdResult_id{1,sti_i};
    curr_edge_select_id = curr_edge_select_id_buff(1:select_neuron_num,:);

    curr_sti_select_neurons_edge_cell = select_neurons_edge_cell{1,sti_i};

    % sort trace by first trial (annotating this means use average sort)
    % first_trial_sti_select_neurons_edge_cell = curr_sti_select_neurons_edge_cell{1,1}(curr_edge_select_id,before+1:end);
    first_trial_sti_select_neurons_edge_cell = curr_sti_select_neurons_edge_cell{1,1}(curr_edge_select_id,:);
    [first_trial_sti_select_neurons_edge_cell_sort,first_sortedTimes, first_sortedValues,first_sortIdx] = sort_trace_by_peakTime(first_trial_sti_select_neurons_edge_cell');
    figure;imagesc(first_trial_sti_select_neurons_edge_cell_sort')
    close all
    figure;
    n =0;
    step = 4;
    % for trial_i = 1:step:size(curr_sti_select_neurons_edge_cell,2)
    for trial_i = 1:5
        n=n+1;
        curr_trial_sti_select_neurons_edge_cell = curr_sti_select_neurons_edge_cell{1,trial_i}(curr_edge_select_id,:);
        curr_trial_sti_select_neurons_edge_cell_sort = curr_trial_sti_select_neurons_edge_cell(first_sortIdx,:);
        % curr_trial_sti_select_neurons_edge_cell_sort = curr_trial_sti_select_neurons_edge_cell(average_sortIdx,:);
        % subplot(1,size(curr_sti_select_neurons_edge_cell,2)/step,n);
        subplot(1,5,n);
        % imagesc(curr_trial_sti_select_neurons_edge_cell_sort)
        imagesc(zscore(curr_trial_sti_select_neurons_edge_cell_sort, 0, 2));hold on;
        colormap(red_blue_color);
        clim([-1,1]);
        y_axis  = 0.5 : 0.5 : size(curr_trial_sti_select_neurons_edge_cell_sort, 1);
        % hold on, plot(2 * ones(1, length(y_axis)), y_axis, '--r', 'linewidth', 2)
        % hold on, plot(before * ones(1, length(y_axis)), y_axis, '--w', 'linewidth', 2)
        % hold on, plot((before+during) * ones(1, length(y_axis)), y_axis, '--r', 'linewidth', 2)
        hold on, plot(12 * ones(12, length(y_axis)), y_axis, '--r', 'linewidth', 2)
        hold on, plot(before_for_corr * ones(1, length(y_axis)), y_axis, '--w', 'linewidth', 2)
        hold on, plot((before_for_corr+during) * ones(1, length(y_axis)), y_axis, '--r', 'linewidth', 2)

        xticks = get(gca, 'XTick'); 
        xticklabels = xticks / param.fs - param.show_sti_before; 
        set(gca, 'XTick', xticks, 'XTickLabel', xticklabels); 

        title(sprintf('trial %d',trial_i), 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold'); 
        xlabel('Time (s)', 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');   
        ylabel('Edge #', 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');    
        set(gca, 'FontSize', 8, 'FontName', 'Arial');            
    end
    c = colorbar('Position', [0.92 0.1 0.02 0.82]); 
    ylabel(c, 'neuron activity intensity value','FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');
    top_num = size(curr_trial_sti_select_neurons_edge_cell_sort,1);
    savefig([figure_path sprintf('multi_session_select_%d_edge_sort_sti_%d_1to5_single_trial',top_num,sti_i)]);
    exportgraphics(gcf,[figure_path sprintf('multi_session_select_%d_edge_sort_sti_%d_1to5_single_trial',top_num,sti_i),'.png'],'Resolution',300)
end

close all
% >>>>>>>>>>>>>>>>>>> sort each single trial all edge
for sti_i = 1:size(select_neurons_edge_cell,2)

    curr_sti_select_neurons_edge_cell = select_neurons_edge_cell{1,sti_i};

    % sort trace by first trial (annotating this means use average sort)
    % first_trial_sti_select_neurons_edge_cell = curr_sti_select_neurons_edge_cell{1,1}(:,before+1:end);
    first_trial_sti_select_neurons_edge_cell = curr_sti_select_neurons_edge_cell{1,1};
    [first_trial_sti_select_neurons_edge_cell_sort,first_sortedTimes, first_sortedValues,first_sortIdx] = sort_trace_by_peakTime(first_trial_sti_select_neurons_edge_cell');
    % figure;imagesc(first_trial_sti_select_neurons_edge_cell_sort')
    % close all
    figure;
    clim_val = [-1 1]; 
    n =0;
    step = 4;
    % for trial_i = 1:step:size(curr_sti_select_neurons_edge_cell,2)
    for trial_i = 1:5
        n=n+1;
        curr_trial_sti_select_neurons_edge_cell = curr_sti_select_neurons_edge_cell{1,trial_i};
        curr_trial_sti_select_neurons_edge_cell_sort = curr_trial_sti_select_neurons_edge_cell(first_sortIdx,:);
        % curr_trial_sti_select_neurons_edge_cell_sort = curr_trial_sti_select_neurons_edge_cell(average_sortIdx,:);
        % subplot(1,size(curr_sti_select_neurons_edge_cell,2)/step,n);
        subplot(1,5,n);
        % imagesc(curr_trial_sti_select_neurons_edge_cell_sort)
        imagesc(zscore(curr_trial_sti_select_neurons_edge_cell_sort, 0, 2));hold on;
        colormap(red_blue_color);
        clim([-1,1]);
        y_axis  = 0.5 : 0.5 : size(curr_trial_sti_select_neurons_edge_cell_sort, 1);
        % hold on, plot(2 * ones(1, length(y_axis)), y_axis, '--r', 'linewidth', 2)
        % hold on, plot(before * ones(1, length(y_axis)), y_axis, '--w', 'linewidth', 2)
        % hold on, plot((before+during) * ones(1, length(y_axis)), y_axis, '--r', 'linewidth', 2)
        hold on, plot(12 * ones(12, length(y_axis)), y_axis, '--r', 'linewidth', 2)
        hold on, plot(before_for_corr * ones(1, length(y_axis)), y_axis, '--w', 'linewidth', 2)
        hold on, plot((before_for_corr+during) * ones(1, length(y_axis)), y_axis, '--r', 'linewidth', 2)

        xticks = get(gca, 'XTick'); 
        xticklabels = xticks / param.fs - param.show_sti_before;
        set(gca, 'XTick', xticks, 'XTickLabel', xticklabels); 

        title(sprintf('trial %d',trial_i), 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold'); 
        xlabel('Time (s)', 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');   
        ylabel('Edge #', 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');    
        set(gca, 'FontSize', 8, 'FontName', 'Arial');            
    end
    c = colorbar('Position', [0.92 0.1 0.02 0.82]); 
    ylabel(c, 'neuron activity intensity value','FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold'); 
    top_num = size(curr_trial_sti_select_neurons_edge_cell_sort,1);
    savefig([figure_path sprintf('multi_session_all_%d_edge_sort_sti_%d_1to5_single_trial',top_num,sti_i)]);
    exportgraphics(gcf,[figure_path sprintf('multi_session_all_%d_edge_sort_sti_%d_1to5_single_trial',top_num,sti_i),'.png'],'Resolution',300)
end



%% check result
color_scheme_npg = get_color();
aa =C_mat(sort_all_neuron_stdResult_id{3}(1:100),:);
[members_trials_cell,stimuli_sort,members_sort_trace,before, during, after, sti_labels,sti_num, trials_num] = sort_multi_neurons_trials_matrix_by_sti(aa, stimuli, param);
deconve_trace_with_label(members_sort_trace, stimuli_sort, color_scheme_npg)

bb =select_neurons_C(sort_select_neuron_stdResult_id{3}(1:47),:);
[members_trials_cell,stimuli_sort,members_sort_trace,before, during, after, sti_labels,sti_num, trials_num] = sort_multi_neurons_trials_matrix_by_sti(bb, stimuli, param);
deconve_trace_with_label(members_sort_trace, stimuli_sort, color_scheme_npg)

cc =select_neurons_C_sort_by_region_corr_matrix(sort_edge_stdResult_id{3}(1:100),:);
[members_trials_cell,stimuli_sort,members_sort_trace,before, during, after, sti_labels,sti_num, trials_num] = sort_multi_neurons_trials_matrix_by_sti(cc, stimuli, param);
deconve_trace_with_label(members_sort_trace, stimuli_sort, color_scheme_npg)



% top TV neurons & ttest neurons & top TV edges
figure('Position', [100, 100, 1500,300]);

set(gcf, 'Color', 'w'); 
numPlots = size(sort_edge_stdResult_value, 2);
for sti_i = 1:numPlots
    subplot(1, numPlots, sti_i); 
try

    % ttest neuron result
    curr_sort_select_neuron_stdResult_value = sort_select_neuron_stdResult_value{sti_i}(1:select_neuron_num);
    [f_neuron, x_neuron] = ksdensity(curr_sort_select_neuron_stdResult_value, 'Support', 'positive'); 
    f_neuron = f_neuron / trapz(x_neuron, f_neuron); 
    h_fill_neuron_ttest = fill([x_neuron, fliplr(x_neuron)], [f_neuron, zeros(size(f_neuron))], [0.2 0.5 0.2], 'FaceAlpha', 0.5, 'EdgeColor', 'none');hold on;
    plot(x_neuron, f_neuron, 'DisplayName', 'Neuron Data', 'LineWidth', 0.1, 'LineStyle', '--','Color', [0.2 0.5 0.2]);hold on;

    % top TV neuron result
    curr_sort_all_neuron_stdResult_value = sort_all_neuron_stdResult_value{sti_i}(1:select_neuron_num);
    [f_neuron, x_neuron] = ksdensity(curr_sort_all_neuron_stdResult_value, 'Support', 'positive'); 
    f_neuron = f_neuron / trapz(x_neuron, f_neuron); 
    h_fill_neuron = fill([x_neuron, fliplr(x_neuron)], [f_neuron, zeros(size(f_neuron))], [0.8 0.8 1], 'FaceAlpha', 0.5, 'EdgeColor', 'none');hold on;
    plot(x_neuron, f_neuron, 'DisplayName', 'Neuron Data', 'LineWidth', 0.8, 'LineStyle', '--','Color', color_scheme_npg(3, :));hold on;

    % top TV edge result
    if process_nozeros == 1
    aa = [];
    aa=nonzeros(sort_edge_stdResult_value{sti_i});
    [f_edge, x_edge] = ksdensity(aa(1:select_neuron_num), 'Support', 'positive'); 
    else
    [f_edge, x_edge] = ksdensity(sort_edge_stdResult_value{sti_i}(1:select_neuron_num), 'Support', 'positive'); 
    end
    f_edge = f_edge / trapz(x_edge, f_edge); 
    h_fill_edge = fill([x_edge, fliplr(x_edge)], [f_edge, zeros(size(f_edge))],  [1 0.85 0.6], 'FaceAlpha', 0.5, 'EdgeColor', 'none');hold on;
    plot(x_edge, f_edge, 'DisplayName', 'Edge Data', 'LineWidth', 0.8,  'Color', color_scheme_npg(2, :)); 


    xlabel('Trial Variance', 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');    
    ylabel('Probability Density', 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');    
    set(gca, 'FontSize', 8, 'FontName', 'Arial');            
    legend([h_fill_neuron_ttest, h_fill_neuron, h_fill_edge], {'Select Neuron','TV Select Neuron', 'TV Select Edge'}, 'Location', 'best', 'FontSize', 6);

    title(sprintf('Stimulus %d', sti_i), 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');
    xlim([0 5]);
    ylim([-0.5,4])
catch
end
end
top_num = size(curr_sort_all_neuron_stdResult_value,1);
savefig(fullfile(figure_path, sprintf('trial_variability_top%dttset_neuron&topTV_neuron&edge_PDF_statistic.fig',top_num)));
exportgraphics(gcf, fullfile(figure_path, sprintf('trial_variability_top%dttset_neuron&topTV_neuron&edge_PDF_statistic.png',top_num)), 'Resolution', 300);


% % >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> average for all sti ,compare sti& blank
% top TV neurons & ttest neurons & top TV edges
figure('Position', [100, 100, 600,300]);


set(gcf, 'Color', 'w'); 
for sti_i = 1:2
    subplot(1, 2, sti_i); 
try

    % ttest neuron result
    curr_sort_select_neuron_stdResult_value = compare_sort_select_neuron_stdResult_value{sti_i}(1:select_neuron_num);
    [f_neuron, x_neuron] = ksdensity(curr_sort_select_neuron_stdResult_value, 'Support', 'positive'); 
    f_neuron = f_neuron / trapz(x_neuron, f_neuron); 
    h_fill_neuron_ttest = fill([x_neuron, fliplr(x_neuron)], [f_neuron, zeros(size(f_neuron))], [0.2 0.5 0.2], 'FaceAlpha', 0.5, 'EdgeColor', 'none');hold on;
    plot(x_neuron, f_neuron, 'DisplayName', 'Neuron Data', 'LineWidth', 0.1, 'LineStyle', '--','Color', [0.2 0.5 0.2]);hold on;

    % top TV neuron result
    curr_sort_all_neuron_stdResult_value = compare_sort_all_neuron_stdResult_value{sti_i}(1:select_neuron_num);
    [f_neuron, x_neuron] = ksdensity(curr_sort_all_neuron_stdResult_value, 'Support', 'positive'); 
    f_neuron = f_neuron / trapz(x_neuron, f_neuron); 
    h_fill_neuron = fill([x_neuron, fliplr(x_neuron)], [f_neuron, zeros(size(f_neuron))], [0.8 0.8 1], 'FaceAlpha', 0.5, 'EdgeColor', 'none');hold on;
    plot(x_neuron, f_neuron, 'DisplayName', 'Neuron Data', 'LineWidth', 0.8, 'LineStyle', '--','Color', color_scheme_npg(3, :));hold on;

    % top TV edge result
    if process_nozeros == 1
    aa = [];
    aa=nonzeros(compare_sort_edge_stdResult_value{sti_i});
    [f_edge, x_edge] = ksdensity(aa(1:select_neuron_num), 'Support', 'positive'); 
    else
    [f_edge, x_edge] = ksdensity(compare_sort_edge_stdResult_value{sti_i}(1:select_neuron_num), 'Support', 'positive');
    end
    f_edge = f_edge / trapz(x_edge, f_edge); 
    h_fill_edge = fill([x_edge, fliplr(x_edge)], [f_edge, zeros(size(f_edge))],  [1 0.85 0.6], 'FaceAlpha', 0.5, 'EdgeColor', 'none');hold on;
    plot(x_edge, f_edge, 'DisplayName', 'Edge Data', 'LineWidth', 0.8,  'Color', color_scheme_npg(2, :)); 

    xlabel('Trial Variance', 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');    
    ylabel('Probability Density', 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');    
    set(gca, 'FontSize', 8, 'FontName', 'Arial');            

    legend([h_fill_neuron_ttest, h_fill_neuron, h_fill_edge], {'Select Neuron','TV Select Neuron', 'TV Select Edge'}, 'Location', 'best', 'FontSize', 6);

    title(sprintf('Stimulus %d', sti_i), 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');
    xlim([0 5]);
    ylim([-0.5,4])
catch
end
end
top_num = size(curr_sort_all_neuron_stdResult_value,1);
savefig(fullfile(figure_path, sprintf('trial_variability_top%dttset_neuron&topTV_neuron&edge_PDF_statistic_average_sti.fig',top_num)));
exportgraphics(gcf, fullfile(figure_path, sprintf('trial_variability_top%dttset_neuron&topTV_neuron&edge_PDF_statistic_average_sti.png',top_num)), 'Resolution', 300);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% subplot
% top TV neurons & ttest neurons & top TV edges
figure('Position', [100, 100, 1500,150]);


set(gcf, 'Color', 'w'); 
numPlots = size(sort_edge_stdResult_value, 2);
for sti_i = 1:numPlots
    subplot(1, numPlots, sti_i); 
    curr_sort_select_neuron_stdResult_value = sort_select_neuron_stdResult_value{sti_i}(1:select_neuron_num);
    [f_neuron, x_neuron] = ksdensity(curr_sort_select_neuron_stdResult_value, 'Support', 'positive'); 
    f_neuron = f_neuron / trapz(x_neuron, f_neuron); 
    h_fill_neuron_ttest = fill([x_neuron, fliplr(x_neuron)], [f_neuron, zeros(size(f_neuron))], [0.2 0.5 0.2], 'FaceAlpha', 0.5, 'EdgeColor', 'none');hold on;
    plot(x_neuron, f_neuron, 'DisplayName', 'Neuron Data', 'LineWidth', 0.8, 'LineStyle', '--','Color', color_scheme_npg(3, :));hold on;

    xlabel('Trial Variance', 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');    
    ylabel('Probability Density', 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');    
    set(gca, 'FontSize', 8, 'FontName', 'Arial');            

    legend([h_fill_neuron_ttest], {'SN'}, 'Location', 'best', 'FontSize', 6);

    title(sprintf('Stimulus %d', sti_i), 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');

end
top_num = size(curr_sort_all_neuron_stdResult_value,1);
savefig(fullfile(figure_path, sprintf('trial_variability_top%dttset_neuron_PDF_statistic.fig',top_num)));
exportgraphics(gcf, fullfile(figure_path, sprintf('trial_variability_top%dttset_neuron_PDF_statistic.png',top_num)), 'Resolution', 300);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% subplot
% top TV neurons & ttest neurons & top TV edges
figure('Position', [100, 100, 600,150]);

set(gcf, 'Color', 'w'); 
numPlots = size(compare_sort_edge_stdResult_value, 2);
for sti_i = 1:numPlots
    subplot(1, numPlots, sti_i); 
    curr_sort_select_neuron_stdResult_value = compare_sort_select_neuron_stdResult_value{sti_i}(1:select_neuron_num);
    [f_neuron, x_neuron] = ksdensity(curr_sort_select_neuron_stdResult_value, 'Support', 'positive'); 
    f_neuron = f_neuron / trapz(x_neuron, f_neuron); 

    h_fill_neuron_ttest = fill([x_neuron, fliplr(x_neuron)], [f_neuron, zeros(size(f_neuron))], [0.2 0.5 0.2], 'FaceAlpha', 0.5, 'EdgeColor', 'none');hold on;
    plot(x_neuron, f_neuron, 'DisplayName', 'Neuron Data', 'LineWidth', 0.8, 'LineStyle', '--','Color', color_scheme_npg(3, :));hold on;

    xlabel('Trial Variance', 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');    
    ylabel('Probability Density', 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');    
    set(gca, 'FontSize', 8, 'FontName', 'Arial');            

    legend([h_fill_neuron_ttest], {'SN'}, 'Location', 'best', 'FontSize', 6);

    title(sprintf('Stimulus %d', sti_i), 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');

end
top_num = size(curr_sort_all_neuron_stdResult_value,1);
savefig(fullfile(figure_path, sprintf('trial_variability_top%dttset_neuron_PDF_statistic_average_sti.fig',top_num)));
exportgraphics(gcf, fullfile(figure_path, sprintf('trial_variability_top%dttset_neuron_PDF_statistic_average_sti.png',top_num)), 'Resolution', 300);



% top TV neurons & ttest neurons & top TV edges
figure('Position', [100, 100, 1500,150]);

set(gcf, 'Color', 'w'); 
numPlots = size(sort_edge_stdResult_value, 2);
for sti_i = 1:numPlots
    subplot(1, numPlots, sti_i); 
try

    if process_nozeros == 1
    aa = [];
    aa = nonzeros(sort_edge_stdResult_value{sti_i});
    [f_edge, x_edge] = ksdensity(aa(1:select_neuron_num), 'Support', 'positive'); 
    else
    [f_edge, x_edge] = ksdensity(sort_edge_stdResult_value{sti_i}(1:select_neuron_num), 'Support', 'positive'); 
    end
    f_edge = f_edge / trapz(x_edge, f_edge); 
    h_fill_edge = fill([x_edge, fliplr(x_edge)], [f_edge, zeros(size(f_edge))],  [1 0.85 0.6], 'FaceAlpha', 0.5, 'EdgeColor', 'none');hold on;
    plot(x_edge, f_edge, 'DisplayName', 'Edge Data', 'LineWidth', 0.8,  'Color', color_scheme_npg(2, :)); 


    xlabel('Trial Variance', 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');   
    ylabel('Probability Density', 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');    
    set(gca, 'FontSize', 8, 'FontName', 'Arial');            

    legend([h_fill_edge], {'SNFC'}, 'Location', 'best', 'FontSize', 6);

    title(sprintf('Stimulus %d', sti_i), 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');

catch
end
end
top_num = size(curr_sort_all_neuron_stdResult_value,1);
savefig(fullfile(figure_path, sprintf('trial_variability_top%dedge_PDF_statistic.fig',top_num)));
exportgraphics(gcf, fullfile(figure_path, sprintf('trial_variability_top%dedge_PDF_statistic.png',top_num)), 'Resolution', 300);

figure('Position', [100, 100, 600,150]);

set(gcf, 'Color', 'w'); 
numPlots = size(compare_sort_edge_stdResult_value, 2);
for sti_i = 1:numPlots
    subplot(1, numPlots, sti_i); 
try

    if process_nozeros == 1
    aa = [];
    aa = nonzeros(compare_sort_edge_stdResult_value{sti_i});
    [f_edge, x_edge] = ksdensity(aa(1:select_neuron_num), 'Support', 'positive'); 
    else
    [f_edge, x_edge] = ksdensity(compare_sort_edge_stdResult_value{sti_i}(1:select_neuron_num), 'Support', 'positive'); 
    end
    f_edge = f_edge / trapz(x_edge, f_edge); 
    h_fill_edge = fill([x_edge, fliplr(x_edge)], [f_edge, zeros(size(f_edge))],  [1 0.85 0.6], 'FaceAlpha', 0.5, 'EdgeColor', 'none');hold on;
    plot(x_edge, f_edge, 'DisplayName', 'Edge Data', 'LineWidth', 0.8,  'Color', color_scheme_npg(2, :)); 


    xlabel('Trial Variance', 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');    
    ylabel('Probability Density', 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');    
    set(gca, 'FontSize', 8, 'FontName', 'Arial');            

    legend([h_fill_edge], {'SNFC'}, 'Location', 'best', 'FontSize', 6);

    title(sprintf('Stimulus %d', sti_i), 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');

catch
end
end
top_num = size(curr_sort_all_neuron_stdResult_value,1);
savefig(fullfile(figure_path, sprintf('trial_variability_top%dedge_PDF_statistic_average_sti.fig',top_num)));
exportgraphics(gcf, fullfile(figure_path, sprintf('trial_variability_top%dedge_PDF_statistic_average_sti.png',top_num)), 'Resolution', 300);





figure('Position', [100, 100, 1500, 350]); 
set(gcf, 'Color', 'w'); 
jitter_amount = 0.15; 
point_size = 36;      
point_alpha = 0.6;    
numPlots = size(sort_select_neuron_stdResult_value, 2);

for sti_i = 1:numPlots
    subplot(1, numPlots, sti_i);
    hold on; 
try
    % 1. T-test neuron result
    curr_ttest_neuron_data = sort_select_neuron_stdResult_value{sti_i}(1:select_neuron_num);
    x_cat1 = ones(size(curr_ttest_neuron_data)); 
    x_jittered1 = x_cat1 + jitter_amount * (rand(size(curr_ttest_neuron_data)) - 0.5) * 2;
    s1 = scatter(x_jittered1, curr_ttest_neuron_data, point_size, ...
                 'filled', 'MarkerFaceColor', [0.2 0.5 0.2], ... % Color from original fill
                 'MarkerFaceAlpha', point_alpha, ...
                 'DisplayName', 'SN');

    % 2. Top TV neuron result
    curr_tv_neuron_data = sort_all_neuron_stdResult_value{sti_i}(1:select_neuron_num);
    x_cat2 = 2 * ones(size(curr_tv_neuron_data)); 
    x_jittered2 = x_cat2 + jitter_amount * (rand(size(curr_tv_neuron_data)) - 0.5) * 2;
    s2 = scatter(x_jittered2, curr_tv_neuron_data, point_size, ...
                 'filled', 'MarkerFaceColor', [0.8 0.8 1], ... % Color from original fill
                 'MarkerFaceAlpha', point_alpha, ...
                 'DisplayName', 'SN(TV)');

    % 3. Top TV edge result
    if process_nozeros == 1
    aa = [];
    aa = nonzeros(sort_edge_stdResult_value{sti_i});
    curr_tv_edge_data = aa(1:select_neuron_num);
    else
    curr_tv_edge_data = sort_edge_stdResult_value{sti_i}(1:select_neuron_num);
    end
    x_cat3 = 3 * ones(size(curr_tv_edge_data)); 
    x_jittered3 = x_cat3 + jitter_amount * (rand(size(curr_tv_edge_data)) - 0.5) * 2;
    s3 = scatter(x_jittered3, curr_tv_edge_data, point_size, ...
                 'filled', 'MarkerFaceColor', [1 0.85 0.6], ... % Color from original fill
                 'MarkerFaceAlpha', point_alpha, ...
                 'DisplayName', 'SNFC(TV)');

    hold off;

    xlabel('Data Group', 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');
    ylabel('Trial Variance', 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold'); 
    set(gca, 'FontSize', 8, 'FontName', 'Arial');

    xticks([1 2 3]);
    xlim([0.5, 3.5]); 

    all_data_current_subplot = [curr_ttest_neuron_data; curr_tv_neuron_data; curr_tv_edge_data];
    min_val = min(all_data_current_subplot);
    max_val = max(all_data_current_subplot);
    padding = (max_val - min_val) * 0.1; % 10% padding
    if min_val == max_val 
        ylim([min_val - 0.1, max_val + 0.1]);
    else
        ylim([min_val - padding, max_val + padding]);
    end
    % ylim([0 2]);

    % 添加图例
    legend([s1, s2, s3], {'SN','SN(TV)', 'SNFC(TV)'}, 'Location', 'best', 'FontSize', 6);

    title(sprintf('Stimulus %d', sti_i), 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');
    grid on; 
    catch
end
end

top_num = select_neuron_num; 
savefig(fullfile(figure_path, sprintf('trial_variability_top%d_ttest_neuron_topTV_neuron_edge_JITTER_SCATTER.fig', top_num)));
exportgraphics(gcf, fullfile(figure_path, sprintf('trial_variability_top%d_ttest_neuron_topTV_neuron_edge_JITTER_SCATTER.png', top_num)), 'Resolution', 300);


figure('Position', [100, 100, 600, 350]); 
set(gcf, 'Color', 'w'); 

jitter_amount = 0.15; 
point_size = 36;      
point_alpha = 0.6;    
numPlots = size(compare_sort_select_neuron_stdResult_value, 2);
for sti_i = 1:numPlots
    subplot(1, numPlots, sti_i);
    hold on; 
try
    % 1. T-test neuron result
    curr_ttest_neuron_data = compare_sort_select_neuron_stdResult_value{sti_i}(1:select_neuron_num);
    x_cat1 = ones(size(curr_ttest_neuron_data)); 
    x_jittered1 = x_cat1 + jitter_amount * (rand(size(curr_ttest_neuron_data)) - 0.5) * 2;
    s1 = scatter(x_jittered1, curr_ttest_neuron_data, point_size, ...
                 'filled', 'MarkerFaceColor', [0.2 0.5 0.2], ... % Color from original fill
                 'MarkerFaceAlpha', point_alpha, ...
                 'DisplayName', 'SN');

    % 2. Top TV neuron result
    curr_tv_neuron_data = compare_sort_all_neuron_stdResult_value{sti_i}(1:select_neuron_num);
    x_cat2 = 2 * ones(size(curr_tv_neuron_data)); 
    x_jittered2 = x_cat2 + jitter_amount * (rand(size(curr_tv_neuron_data)) - 0.5) * 2;
    s2 = scatter(x_jittered2, curr_tv_neuron_data, point_size, ...
                 'filled', 'MarkerFaceColor', [0.8 0.8 1], ... % Color from original fill
                 'MarkerFaceAlpha', point_alpha, ...
                 'DisplayName', 'SN(TV)');

    % 3. Top TV edge result
    if process_nozeros == 1
    aa = [];
    aa = nonzeros(compare_sort_edge_stdResult_value{sti_i});
    curr_tv_edge_data = aa(1:select_neuron_num);
    else
    curr_tv_edge_data = compare_sort_edge_stdResult_value{sti_i}(1:select_neuron_num);
    end
    x_cat3 = 3 * ones(size(curr_tv_edge_data)); 
    x_jittered3 = x_cat3 + jitter_amount * (rand(size(curr_tv_edge_data)) - 0.5) * 2;
    s3 = scatter(x_jittered3, curr_tv_edge_data, point_size, ...
                 'filled', 'MarkerFaceColor', [1 0.85 0.6], ... % Color from original fill
                 'MarkerFaceAlpha', point_alpha, ...
                 'DisplayName', 'SNFC(TV)');

    hold off;

    xlabel('Data Group', 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');
    ylabel('Trial Variance', 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold'); % Y轴标签更新
    set(gca, 'FontSize', 8, 'FontName', 'Arial');


    xticks([1 2 3]);
    xlim([0.5, 3.5]); 

    all_data_current_subplot = [curr_ttest_neuron_data; curr_tv_neuron_data; curr_tv_edge_data];
    min_val = min(all_data_current_subplot);
    max_val = max(all_data_current_subplot);
    padding = (max_val - min_val) * 0.1; % 10% padding
    if min_val == max_val 
        ylim([min_val - 0.1, max_val + 0.1]);
    else
        ylim([min_val - padding, max_val + padding]);
    end
    legend([s1, s2, s3], {'SN','SN(TV)', 'SNFC(TV)'}, 'Location', 'best', 'FontSize', 6);

    title(sprintf('Stimulus %d', sti_i), 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');
    grid on; 
    catch
end
end

top_num = select_neuron_num; 
savefig(fullfile(figure_path, sprintf('trial_variability_top%d_ttest_neuron_topTV_neuron_edge_JITTER_SCATTER_average_sti.fig', top_num)));
exportgraphics(gcf, fullfile(figure_path, sprintf('trial_variability_top%d_ttest_neuron_topTV_neuron_edge_JITTER_SCATTER_average_sti.png', top_num)), 'Resolution', 300);



figure('Position', [100, 100, 1500, 350]); 
set(gcf, 'Color', 'w'); 

jitter_amount = 0.15; 
point_size = 36;     
point_alpha = 0.6;    
numPlots = size(sort_select_neuron_stdResult_value, 2);
for sti_i = 1:numPlots
    subplot(1, numPlots, sti_i);
    hold on; 
try
    % 1. T-test neuron result
    curr_ttest_neuron_data = sort_select_neuron_stdResult_value{sti_i}(1:select_neuron_num);
    x_cat1 = ones(size(curr_ttest_neuron_data)); 
    x_jittered1 = x_cat1 + jitter_amount * (rand(size(curr_ttest_neuron_data)) - 0.5) * 2;
    s1 = scatter(x_jittered1, curr_ttest_neuron_data, point_size, ...
                 'filled', 'MarkerFaceColor', [0.2 0.5 0.2], ... % Color from original fill
                 'MarkerFaceAlpha', point_alpha, ...
                 'DisplayName', 'SN');

    % % 2. Top TV neuron result
    % curr_tv_neuron_data = sort_all_neuron_stdResult_value{sti_i}(1:select_neuron_num);
    % x_cat2 = 2 * ones(size(curr_tv_neuron_data)); 
    % x_jittered2 = x_cat2 + jitter_amount * (rand(size(curr_tv_neuron_data)) - 0.5) * 2;
    % s2 = scatter(x_jittered2, curr_tv_neuron_data, point_size, ...
    %              'filled', 'MarkerFaceColor', [0.8 0.8 1], ... % Color from original fill
    %              'MarkerFaceAlpha', point_alpha, ...
    %              'DisplayName', 'TV Neuron');

    % 3. Top TV edge result
    curr_tv_edge_data = sort_edge_stdResult_value{sti_i}(1:select_neuron_num);
    x_cat3 = 2 * ones(size(curr_tv_edge_data)); 
    x_jittered3 = x_cat3 + jitter_amount * (rand(size(curr_tv_edge_data)) - 0.5) * 2;
    s3 = scatter(x_jittered3, curr_tv_edge_data, point_size, ...
                 'filled', 'MarkerFaceColor', [1 0.85 0.6], ... % Color from original fill
                 'MarkerFaceAlpha', point_alpha, ...
                 'DisplayName', 'SNFC');

    hold off;
    xlabel('Data Group', 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');
    ylabel('Trial Variance', 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold'); 
    set(gca, 'FontSize', 8, 'FontName', 'Arial');
    xticks([1 2]);
    xlim([0.5, 2.5]);
    all_data_current_subplot = [curr_ttest_neuron_data; curr_tv_neuron_data; curr_tv_edge_data];
    min_val = min(all_data_current_subplot);
    max_val = max(all_data_current_subplot);
    padding = (max_val - min_val) * 0.1; % 10% padding
    if min_val == max_val 
        ylim([min_val - 0.1, max_val + 0.1]);
    else
        ylim([min_val - padding, max_val + padding]);
    end
    legend([s1, s2, s3], {'SN', 'SNFC'}, 'Location', 'best', 'FontSize', 6);

    title(sprintf('Stimulus %d', sti_i), 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');
    grid on; 
    catch
end
end

top_num = select_neuron_num; 
savefig(fullfile(figure_path, sprintf('top%d_ttest_neuron_topTV_edge_JITTER_SCATTER.fig', top_num)));
exportgraphics(gcf, fullfile(figure_path, sprintf('top%d_ttest_neuron_topTV_edge_JITTER_SCATTER.png', top_num)), 'Resolution', 300);


figure('Position', [100, 100, 600, 350]);
set(gcf, 'Color', 'w'); 
jitter_amount = 0.15; 
point_size = 36;      
point_alpha = 0.6;    
numPlots = size(compare_sort_select_neuron_stdResult_value, 2);
for sti_i = 1:numPlots
    subplot(1, numPlots, sti_i);
    hold on; 
try
    % 1. T-test neuron result
    curr_ttest_neuron_data = compare_sort_select_neuron_stdResult_value{sti_i}(1:select_neuron_num);
    x_cat1 = ones(size(curr_ttest_neuron_data)); 
    x_jittered1 = x_cat1 + jitter_amount * (rand(size(curr_ttest_neuron_data)) - 0.5) * 2;
    s1 = scatter(x_jittered1, curr_ttest_neuron_data, point_size, ...
                 'filled', 'MarkerFaceColor', [0.2 0.5 0.2], ... % Color from original fill
                 'MarkerFaceAlpha', point_alpha, ...
                 'DisplayName', 'SN');

    % % 2. Top TV neuron result
    % curr_tv_neuron_data = sort_all_neuron_stdResult_value{sti_i}(1:select_neuron_num);
    % x_cat2 = 2 * ones(size(curr_tv_neuron_data)); 
    % x_jittered2 = x_cat2 + jitter_amount * (rand(size(curr_tv_neuron_data)) - 0.5) * 2;
    % s2 = scatter(x_jittered2, curr_tv_neuron_data, point_size, ...
    %              'filled', 'MarkerFaceColor', [0.8 0.8 1], ... % Color from original fill
    %              'MarkerFaceAlpha', point_alpha, ...
    %              'DisplayName', 'TV Neuron');

    % 3. Top TV edge result
    curr_tv_edge_data = compare_sort_edge_stdResult_value{sti_i}(1:select_neuron_num);
    x_cat3 = 2 * ones(size(curr_tv_edge_data)); 
    x_jittered3 = x_cat3 + jitter_amount * (rand(size(curr_tv_edge_data)) - 0.5) * 2;
    s3 = scatter(x_jittered3, curr_tv_edge_data, point_size, ...
                 'filled', 'MarkerFaceColor', [1 0.85 0.6], ... % Color from original fill
                 'MarkerFaceAlpha', point_alpha, ...
                 'DisplayName', 'SNFC');

    hold off;
    xlabel('Data Group', 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');
    ylabel('Trial Variance', 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold'); % Y轴标签更新
    set(gca, 'FontSize', 8, 'FontName', 'Arial');

    xticks([1 2]);
    xlim([0.5, 2.5]); 

    all_data_current_subplot = [curr_ttest_neuron_data; curr_tv_neuron_data; curr_tv_edge_data];
    min_val = min(all_data_current_subplot);
    max_val = max(all_data_current_subplot);
    padding = (max_val - min_val) * 0.1; % 10% padding
    if min_val == max_val 
        ylim([min_val - 0.1, max_val + 0.1]);
    else
        ylim([min_val - padding, max_val + padding]);
    end
    legend([s1, s2, s3], {'SN', 'SNFC'}, 'Location', 'best', 'FontSize', 6);

    title(sprintf('Stimulus %d', sti_i), 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');
    grid on;
    catch
end
end

top_num = select_neuron_num;
savefig(fullfile(figure_path, sprintf('top%d_ttest_neuron_topTV_edge_JITTER_SCATTER_average_sti.fig', top_num)));
exportgraphics(gcf, fullfile(figure_path, sprintf('top%d_ttest_neuron_topTV_edge_JITTER_SCATTER_average_sti.png', top_num)), 'Resolution', 300);



% top TV neurons & ttest neurons & top TV edges
figure('Position', [100, 100, 1500, 300]);
set(gcf, 'Color', 'w'); 
numPlots = size(sort_edge_stdResult_value, 2);
for sti_i = 1:numPlots
    subplot(1, numPlots, sti_i);
try
    % ttest neuron result
    curr_sort_select_neuron_stdResult_value = sort_select_neuron_stdResult_value{sti_i}(1:select_neuron_num);
    [f_neuron_ttest, x_neuron_ttest] = ksdensity(curr_sort_select_neuron_stdResult_value, 'Support', 'positive'); 
    f_neuron_ttest = f_neuron_ttest / trapz(x_neuron_ttest, f_neuron_ttest); 
    h_fill_neuron_ttest = fill([x_neuron_ttest, fliplr(x_neuron_ttest)], [f_neuron_ttest, zeros(size(f_neuron_ttest))], [0.2 0.5 0.2], 'FaceAlpha', 0.5, 'EdgeColor', 'none');hold on;
    plot(x_neuron_ttest, f_neuron_ttest, 'DisplayName', 'SN', 'LineWidth', 0.8, 'LineStyle', '--','Color', color_scheme_npg(3, :));hold on; 

    % top TV neuron result
    % 计算并绘制 sort_select_neuron_stdResult_value 的核密度估计
    curr_sort_all_neuron_stdResult_value = sort_all_neuron_stdResult_value{sti_i}(1:select_neuron_num);
    [f_neuron_tv, x_neuron_tv] = ksdensity(curr_sort_all_neuron_stdResult_value, 'Support', 'positive'); 
    f_neuron_tv = f_neuron_tv / trapz(x_neuron_tv, f_neuron_tv); 
    h_fill_neuron = fill([x_neuron_tv, fliplr(x_neuron_tv)], [f_neuron_tv, zeros(size(f_neuron_tv))], [0.8 0.8 1], 'FaceAlpha', 0.5, 'EdgeColor', 'none');hold on;
    plot(x_neuron_tv, f_neuron_tv, 'DisplayName', 'SN(TV)', 'LineWidth', 0.8, 'LineStyle', '--','Color', color_scheme_npg(3, :));hold on; 

    % top TV edge result
    if process_nozeros == 1
    aa = [];
    aa = nonzeros(sort_edge_stdResult_value{sti_i});
    [f_edge, x_edge] = ksdensity(aa(1:select_neuron_num), 'Support', 'positive');
    else
    [f_edge, x_edge] = ksdensity(sort_edge_stdResult_value{sti_i}(1:select_neuron_num), 'Support', 'positive'); % 核密度估计
    end
    f_edge = f_edge / trapz(x_edge, f_edge); % 确保密度总和为1
    h_fill_edge = fill([x_edge, fliplr(x_edge)], [f_edge, zeros(size(f_edge))],  [1 0.85 0.6], 'FaceAlpha', 0.5, 'EdgeColor', 'none');hold on;
    plot(x_edge, f_edge, 'DisplayName', 'SNFC (TV)', 'LineWidth', 0.8,  'Color', color_scheme_npg(2, :)); 
    set(gca, 'XScale', 'log'); % 将X轴设置为对数刻度

    xlabel('Trial Variance (log scale)', 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');    
    ylabel('Probability Density', 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');    
    set(gca, 'FontSize', 8, 'FontName', 'Arial');           
    legend([h_fill_neuron_ttest, h_fill_neuron, h_fill_edge], {'SN','SN(TV)', 'SNFC'}, 'Location', 'best', 'FontSize', 6);

    title(sprintf('Stimulus %d', sti_i), 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');
    all_x_current_plot = [x_neuron_ttest(x_neuron_ttest>0), x_neuron_tv(x_neuron_tv>0), x_edge(x_edge>0)];
    if ~isempty(all_x_current_plot)
        min_x_val = min(all_x_current_plot);
        dynamic_min_xlim = max(min_x_val * 0.5, 1e-4); 
        xlim([dynamic_min_xlim 1000]);
    else
        xlim([0.001 2]); % 如果没有有效数据点，使用一个默认的小正数
    end
    ylim([-0.5,10])
    if sti_i == 5
        xlim([dynamic_min_xlim 10000]);
        ylim([-0.2,3])
    end 
    catch
end
end
top_num = size(curr_sort_all_neuron_stdResult_value,1); 
savefig(fullfile(figure_path, sprintf('trial_variability_top%dttset_neuron&topTV_neuron&edge_PDF_statistic_logX.fig',select_neuron_num))); 
exportgraphics(gcf, fullfile(figure_path, sprintf('trial_variability_top%dttset_neuron&topTV_neuron&edge_PDF_statistic_logX.png',select_neuron_num)), 'Resolution', 300); 


figure('Position', [100, 100, 600, 300]);
set(gcf, 'Color', 'w'); 
numPlots = size(compare_sort_edge_stdResult_value, 2);
for sti_i = 1:numPlots
    subplot(1, numPlots, sti_i);
try
    % ttest neuron result
    curr_sort_select_neuron_stdResult_value = compare_sort_select_neuron_stdResult_value{sti_i}(1:select_neuron_num);
    [f_neuron_ttest, x_neuron_ttest] = ksdensity(curr_sort_select_neuron_stdResult_value, 'Support', 'positive'); 
    f_neuron_ttest = f_neuron_ttest / trapz(x_neuron_ttest, f_neuron_ttest); 
    h_fill_neuron_ttest = fill([x_neuron_ttest, fliplr(x_neuron_ttest)], [f_neuron_ttest, zeros(size(f_neuron_ttest))], [0.2 0.5 0.2], 'FaceAlpha', 0.5, 'EdgeColor', 'none');hold on;
    plot(x_neuron_ttest, f_neuron_ttest, 'DisplayName', 'SN', 'LineWidth', 0.8, 'LineStyle', '--','Color', color_scheme_npg(3, :));hold on; 

    % top TV neuron result
    % 计算并绘制 sort_select_neuron_stdResult_value 的核密度估计
    curr_sort_all_neuron_stdResult_value = compare_sort_all_neuron_stdResult_value{sti_i}(1:select_neuron_num);
    [f_neuron_tv, x_neuron_tv] = ksdensity(curr_sort_all_neuron_stdResult_value, 'Support', 'positive'); 
    f_neuron_tv = f_neuron_tv / trapz(x_neuron_tv, f_neuron_tv);
    h_fill_neuron = fill([x_neuron_tv, fliplr(x_neuron_tv)], [f_neuron_tv, zeros(size(f_neuron_tv))], [0.8 0.8 1], 'FaceAlpha', 0.5, 'EdgeColor', 'none');hold on;
    plot(x_neuron_tv, f_neuron_tv, 'DisplayName', 'SN(TV)', 'LineWidth', 0.8, 'LineStyle', '--','Color', color_scheme_npg(3, :));hold on; 

    % top TV edge result
    % 计算并绘制 sort_edge_stdResult_value 的核密度估计
    if process_nozeros == 1
    aa = [];
    aa = nonzeros(compare_sort_edge_stdResult_value{sti_i});
    [f_edge, x_edge] = ksdensity(aa(1:select_neuron_num), 'Support', 'positive');
    else
    [f_edge, x_edge] = ksdensity(compare_sort_edge_stdResult_value{sti_i}(1:select_neuron_num), 'Support', 'positive'); 
    end
    f_edge = f_edge / trapz(x_edge, f_edge); % 确保密度总和为1
    h_fill_edge = fill([x_edge, fliplr(x_edge)], [f_edge, zeros(size(f_edge))],  [1 0.85 0.6], 'FaceAlpha', 0.5, 'EdgeColor', 'none');hold on;
    plot(x_edge, f_edge, 'DisplayName', 'SNFC (TV)', 'LineWidth', 0.8,  'Color', color_scheme_npg(2, :)); 

    set(gca, 'XScale', 'log'); % 将X轴设置为对数刻度
    xlabel('Trial Variance (log scale)', 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');    
    ylabel('Probability Density', 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');    
    set(gca, 'FontSize', 8, 'FontName', 'Arial');            
    legend([h_fill_neuron_ttest, h_fill_neuron, h_fill_edge], {'SN','SN(TV)', 'SNFC'}, 'Location', 'best', 'FontSize', 6);

    title(sprintf('Stimulus %d', sti_i), 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');
    all_x_current_plot = [x_neuron_ttest(x_neuron_ttest>0), x_neuron_tv(x_neuron_tv>0), x_edge(x_edge>0)];
    if ~isempty(all_x_current_plot)
        min_x_val = min(all_x_current_plot);
        dynamic_min_xlim = max(min_x_val * 0.5, 1e-4); % 避免0或负数，并设置一个绝对最小界限
        xlim([dynamic_min_xlim 1000]);
    else
        % xlim([0.001 2]); 
    end
    ylim([-0.5,10])
    if sti_i == 5
        xlim([dynamic_min_xlim 10000]);
        ylim([-0.2,3])
    end 
    catch
end
end
top_num = size(curr_sort_all_neuron_stdResult_value,1); 
savefig(fullfile(figure_path, sprintf('trial_variability_top%dttset_neuron&topTV_neuron&edge_PDF_statistic_logX_average_sti.fig',select_neuron_num))); 
exportgraphics(gcf, fullfile(figure_path, sprintf('trial_variability_top%dttset_neuron&topTV_neuron&edge_PDF_statistic_logX_average_sti.png',select_neuron_num)), 'Resolution', 300); % 修改文件名



%% PDF plot

figure('Position', [100, 100, 1500,150]);
set(gcf, 'Color', 'w'); 
numPlots = size(sort_edge_stdResult_value, 2);
for sti_i = 1:numPlots
    subplot(1, numPlots, sti_i); 
    try

    % 计算并绘制 sort_all_neuron_stdResult_value 的核密度估计
    % [f_neuron, x_neuron]= ksdensity(sort_all_neuron_stdResult_value{sti_i}(1:100), 'Function', 'cdf'); % 累计分布函数 (CDF)
    curr_sort_all_neuron_stdResult_value = sort_all_neuron_stdResult_value{sti_i}(1:select_neuron_num);
    [f_neuron, x_neuron] = ksdensity(curr_sort_all_neuron_stdResult_value, 'Support', 'positive'); 
    f_neuron = f_neuron / trapz(x_neuron, f_neuron); 
    h_fill_neuron = fill([x_neuron, fliplr(x_neuron)], [f_neuron, zeros(size(f_neuron))], [0.8 0.8 1], 'FaceAlpha', 0.5, 'EdgeColor', 'none');hold on;
    plot(x_neuron, f_neuron, 'DisplayName', 'Neuron Data', 'LineWidth', 0.8, 'LineStyle', '--','Color', color_scheme_npg(3, :));hold on;

    % 计算并绘制 sort_edge_stdResult_value 的核密度估计
    if process_nozeros == 1
    aa = [];
    aa = nonzeros(sort_edge_stdResult_value{sti_i});
    [f_edge, x_edge] = ksdensity(aa(1:select_neuron_num), 'Support', 'positive');
    else
    [f_edge, x_edge] = ksdensity(sort_edge_stdResult_value{sti_i}(1:select_neuron_num), 'Support', 'positive'); 
    end
    f_edge = f_edge / trapz(x_edge, f_edge); 
    h_fill_edge = fill([x_edge, fliplr(x_edge)], [f_edge, zeros(size(f_edge))],  [1 0.85 0.6], 'FaceAlpha', 0.5, 'EdgeColor', 'none');hold on;
    plot(x_edge, f_edge, 'DisplayName', 'Edge Data', 'LineWidth', 0.8,  'Color', color_scheme_npg(2, :)); 

    xlabel('Trial Variance', 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');    
    ylabel('Probability Density', 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');   
    set(gca, 'FontSize', 8, 'FontName', 'Arial');           
    legend([h_fill_neuron, h_fill_edge], {'SN(TV)', 'SNFC'}, 'Location', 'best', 'FontSize', 6);

    title(sprintf('Stimulus %d', sti_i), 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');
    catch
end

end
top_num = size(curr_sort_all_neuron_stdResult_value,1);
savefig(fullfile(figure_path, sprintf('trial_variability_top%d_neuron&edge_PDF_statistic.fig',top_num)));
exportgraphics(gcf, fullfile(figure_path, sprintf('trial_variability_top%d_neuron&edge_PDF_statistic.png',top_num)), 'Resolution', 300); % todo

figure('Position', [100, 100, 600,150]);
set(gcf, 'Color', 'w'); 
numPlots = size(compare_sort_edge_stdResult_value, 2);
for sti_i = 1:numPlots
    subplot(1, numPlots, sti_i); 
    try

    % 计算并绘制 sort_all_neuron_stdResult_value 的核密度估计
    % [f_neuron, x_neuron]= ksdensity(sort_all_neuron_stdResult_value{sti_i}(1:100), 'Function', 'cdf'); % 累计分布函数 (CDF)
    curr_sort_all_neuron_stdResult_value = compare_sort_all_neuron_stdResult_value{sti_i}(1:select_neuron_num);
    [f_neuron, x_neuron] = ksdensity(curr_sort_all_neuron_stdResult_value, 'Support', 'positive'); % 核密度估计 % 概率密度函数 (PDF)
    f_neuron = f_neuron / trapz(x_neuron, f_neuron); 
    h_fill_neuron = fill([x_neuron, fliplr(x_neuron)], [f_neuron, zeros(size(f_neuron))], [0.8 0.8 1], 'FaceAlpha', 0.5, 'EdgeColor', 'none');hold on;
    plot(x_neuron, f_neuron, 'DisplayName', 'Neuron Data', 'LineWidth', 0.8, 'LineStyle', '--','Color', color_scheme_npg(3, :));hold on;

    if process_nozeros == 1
    aa = [];
    aa = nonzeros(compare_sort_edge_stdResult_value{sti_i});
    [f_edge, x_edge] = ksdensity(aa(1:select_neuron_num), 'Support', 'positive');
    else
    [f_edge, x_edge] = ksdensity(compare_sort_edge_stdResult_value{sti_i}(1:select_neuron_num), 'Support', 'positive'); 
    end
    f_edge = f_edge / trapz(x_edge, f_edge); 
    h_fill_edge = fill([x_edge, fliplr(x_edge)], [f_edge, zeros(size(f_edge))],  [1 0.85 0.6], 'FaceAlpha', 0.5, 'EdgeColor', 'none');hold on;
    plot(x_edge, f_edge, 'DisplayName', 'Edge Data', 'LineWidth', 0.8,  'Color', color_scheme_npg(2, :)); 


    xlabel('Trial Variance', 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');    
    ylabel('Probability Density', 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');    
    set(gca, 'FontSize', 8, 'FontName', 'Arial');            

    legend([h_fill_neuron, h_fill_edge], {'SN(TV)', 'SNFC'}, 'Location', 'best', 'FontSize', 6);

    title(sprintf('Stimulus %d', sti_i), 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');
    catch
end

end
top_num = size(curr_sort_all_neuron_stdResult_value,1);
savefig(fullfile(figure_path, sprintf('trial_variability_top%d_neuron&edge_PDF_statistic_average_sti.fig',top_num)));
exportgraphics(gcf, fullfile(figure_path, sprintf('trial_variability_top%d_neuron&edge_PDF_statistic_average_sti.png',top_num)), 'Resolution', 300); % todo


% % select neuron & edge top TVV
figure('Position', [100, 100, 1500,150]);
for sti_i = 1:size(sort_edge_stdResult_value,2)
    subplot(1,size(sort_edge_stdResult_value,2),sti_i)
try
    curr_sort_select_neuron_stdResult_value = sort_select_neuron_stdResult_value{sti_i}(1:select_neuron_num);
    [f_neuron, x_neuron] = ksdensity(curr_sort_select_neuron_stdResult_value, 'Support', 'positive'); 
    f_neuron = f_neuron / trapz(x_neuron, f_neuron); 
    h_fill_neuron = fill([x_neuron, fliplr(x_neuron)], [f_neuron, zeros(size(f_neuron))], [0.8 0.8 1], 'FaceAlpha', 0.5, 'EdgeColor', 'none');hold on;
    plot(x_neuron, f_neuron, 'DisplayName', 'Neuron Data', 'LineWidth', 0.8, 'LineStyle', '--','Color', color_scheme_npg(3, :));hold on;

    if process_nozeros == 1
    aa = [];
    aa = nonzeros(sort_edge_stdResult_value{sti_i});
    [f_edge, x_edge] = ksdensity(aa(1:select_neuron_num), 'Support', 'positive');
    else
    [f_edge, x_edge] = ksdensity(sort_edge_stdResult_value{sti_i}(1:select_neuron_num), 'Support', 'positive'); 
    end
    f_edge = f_edge / trapz(x_edge, f_edge); 
    h_fill_edge = fill([x_edge, fliplr(x_edge)], [f_edge, zeros(size(f_edge))],  [1 0.85 0.6], 'FaceAlpha', 0.5, 'EdgeColor', 'none');hold on;
    plot(x_edge, f_edge, 'DisplayName', 'Edge Data', 'LineWidth', 0.8,  'Color', color_scheme_npg(2, :)); 

    xlim([0 10]);

    xlabel('Trial Variance', 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');    
    ylabel('Probability Density', 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');    
    set(gca, 'FontSize', 8, 'FontName', 'Arial');           

    legend([h_fill_neuron, h_fill_edge], {'SN', 'SNFC'}, 'Location', 'best', 'FontSize', 6);
    title(sprintf('Stimulus %d', sti_i), 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');
catch
end
end
top_num = size(curr_sort_select_neuron_stdResult_value,1);

savefig(fullfile(figure_path, sprintf('trial_variability_ttest_top%d_neuron&edge_PDF_statistic.fig',top_num)));
exportgraphics(gcf, fullfile(figure_path, sprintf('trial_variability_ttest_top%d_neuron&edge_PDF_statistic.png',top_num)), 'Resolution', 300);


% % select neuron & edge top TVV
figure('Position', [100, 100, 600,150]);
for sti_i = 1:size(compare_sort_edge_stdResult_value,2)
    subplot(1,size(compare_sort_edge_stdResult_value,2),sti_i)
try
    % 计算密度估计
    curr_sort_select_neuron_stdResult_value = compare_sort_select_neuron_stdResult_value{sti_i}(1:select_neuron_num);
    [f_neuron, x_neuron] = ksdensity(curr_sort_select_neuron_stdResult_value, 'Support', 'positive'); 
    f_neuron = f_neuron / trapz(x_neuron, f_neuron); 
    h_fill_neuron = fill([x_neuron, fliplr(x_neuron)], [f_neuron, zeros(size(f_neuron))], [0.8 0.8 1], 'FaceAlpha', 0.5, 'EdgeColor', 'none');hold on;
    plot(x_neuron, f_neuron, 'DisplayName', 'Neuron Data', 'LineWidth', 0.8, 'LineStyle', '--','Color', color_scheme_npg(3, :));hold on;

    % 计算密度估计
    if process_nozeros == 1
    aa = [];
    aa = nonzeros(compare_sort_edge_stdResult_value{sti_i});
    [f_edge, x_edge] = ksdensity(aa(1:select_neuron_num), 'Support', 'positive');
    else
    [f_edge, x_edge] = ksdensity(compare_sort_edge_stdResult_value{sti_i}(1:select_neuron_num), 'Support', 'positive'); 
    end
    f_edge = f_edge / trapz(x_edge, f_edge); % 确保密度总和为1
    h_fill_edge = fill([x_edge, fliplr(x_edge)], [f_edge, zeros(size(f_edge))],  [1 0.85 0.6], 'FaceAlpha', 0.5, 'EdgeColor', 'none');hold on;
    plot(x_edge, f_edge, 'DisplayName', 'Edge Data', 'LineWidth', 0.8,  'Color', color_scheme_npg(2, :)); 

    xlim([0 10]);

    xlabel('Trial Variance', 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');    
    ylabel('Probability Density', 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');    
    set(gca, 'FontSize', 8, 'FontName', 'Arial');           

    legend([h_fill_neuron, h_fill_edge], {'SN', 'SNFC'}, 'Location', 'best', 'FontSize', 6);
    title(sprintf('Stimulus %d', sti_i), 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');
catch
end
end
top_num = size(curr_sort_select_neuron_stdResult_value,1);

savefig(fullfile(figure_path, sprintf('trial_variability_ttest_top%d_neuron&edge_PDF_statistic_average_sti.fig',top_num)));
exportgraphics(gcf, fullfile(figure_path, sprintf('trial_variability_ttest_top%d_neuron&edge_PDF_statistic_average_sti.png',top_num)), 'Resolution', 300);


% % select neuron & all edge top TVV
figure('Position', [100, 100, 1500,150]);

for sti_i = 1:size(sort_edge_stdResult_value,2)
    subplot(1,size(sort_edge_stdResult_value,2),sti_i)
try
    % 计算密度估计
    curr_sort_all_neuron_stdResult_value = sort_all_neuron_stdResult_value{sti_i};
    [f_neuron, x_neuron] = ksdensity(curr_sort_all_neuron_stdResult_value, 'Support', 'positive'); % 使用核密度估计
    f_neuron = f_neuron / trapz(x_neuron, f_neuron); % 确保密度总和为1
    h_fill_neuron = fill([x_neuron, fliplr(x_neuron)], [f_neuron, zeros(size(f_neuron))], [0.8 0.8 1], 'FaceAlpha', 0.5, 'EdgeColor', 'none');hold on;
    plot(x_neuron, f_neuron, 'DisplayName', 'Neuron Data', 'LineWidth', 0.8, 'LineStyle', '--','Color', color_scheme_npg(3, :));hold on;

    % 计算密度估计
    if process_nozeros == 1
    aa = [];
    aa = nonzeros(sort_edge_stdResult_value{sti_i});
    % [f_edge, x_edge] = ksdensity(aa(1:select_neuron_num), 'Support', 'positive');
    [f_edge, x_edge] = ksdensity(aa, 'Support', 'positive');
    else
    [f_edge, x_edge] = ksdensity(sort_edge_stdResult_value{sti_i}, 'Support', 'positive'); 
    end
    f_edge = f_edge / trapz(x_edge, f_edge); 
    h_fill_edge = fill([x_edge, fliplr(x_edge)], [f_edge, zeros(size(f_edge))],  [1 0.85 0.6], 'FaceAlpha', 0.5, 'EdgeColor', 'none');hold on;
    plot(x_edge, f_edge, 'DisplayName', 'Edge Data', 'LineWidth', 0.8,  'Color', color_scheme_npg(2, :)); 
    xlim([-10 200]);
    xlabel('Trial Variance', 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');    
    ylabel('Probability Density', 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');    
    set(gca, 'FontSize', 8, 'FontName', 'Arial');            

    legend([h_fill_neuron, h_fill_edge], {'All Neuron', 'All Edge'}, 'Location', 'best', 'FontSize', 6);
    title(sprintf('Stimulus %d', sti_i), 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');
catch
end
end

savefig([figure_path sprintf('trial_variability_all_neuron&edge_PDF_statistic')]);
exportgraphics(gcf,[figure_path sprintf('trial_variability_all_neuron&edge_PDF_statistic'),'.png'],'Resolution',300)



figure('Position', [100, 100, 600,150]);

for sti_i = 1:size(compare_sort_edge_stdResult_value,2)
    subplot(1,size(compare_sort_edge_stdResult_value,2),sti_i)
try
    % 计算密度估计
    curr_sort_all_neuron_stdResult_value = compare_sort_all_neuron_stdResult_value{sti_i};
    [f_neuron, x_neuron] = ksdensity(curr_sort_all_neuron_stdResult_value, 'Support', 'positive'); 
    f_neuron = f_neuron / trapz(x_neuron, f_neuron); 
    h_fill_neuron = fill([x_neuron, fliplr(x_neuron)], [f_neuron, zeros(size(f_neuron))], [0.8 0.8 1], 'FaceAlpha', 0.5, 'EdgeColor', 'none');hold on;
    plot(x_neuron, f_neuron, 'DisplayName', 'Neuron Data', 'LineWidth', 0.8, 'LineStyle', '--','Color', color_scheme_npg(3, :));hold on;

    if process_nozeros == 1
    aa = [];
    aa = nonzeros(compare_sort_edge_stdResult_value{sti_i});
    % [f_edge, x_edge] = ksdensity(aa(1:select_neuron_num), 'Support', 'positive');
    [f_edge, x_edge] = ksdensity(aa, 'Support', 'positive');
    else
    [f_edge, x_edge] = ksdensity(compare_sort_edge_stdResult_value{sti_i}, 'Support', 'positive'); 
    end
    f_edge = f_edge / trapz(x_edge, f_edge); 
    h_fill_edge = fill([x_edge, fliplr(x_edge)], [f_edge, zeros(size(f_edge))],  [1 0.85 0.6], 'FaceAlpha', 0.5, 'EdgeColor', 'none');hold on;
    plot(x_edge, f_edge, 'DisplayName', 'Edge Data', 'LineWidth', 0.8,  'Color', color_scheme_npg(2, :)); 
    xlim([-10 200]);
    xlabel('Trial Variance', 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');    
    ylabel('Probability Density', 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');    
    set(gca, 'FontSize', 8, 'FontName', 'Arial');            

    % 添加图例与样式
    legend([h_fill_neuron, h_fill_edge], {'All Neuron', 'All Edge'}, 'Location', 'best', 'FontSize', 6);
    title(sprintf('Stimulus %d', sti_i), 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');
catch
end
end

savefig([figure_path sprintf('trial_variability_all_neuron&edge_PDF_statistic_average_sti')]);
exportgraphics(gcf,[figure_path sprintf('trial_variability_all_neuron&edge_PDF_statistic_average_sti'),'.png'],'Resolution',300)





% % select neuron & all edge top TVV
figure('Position', [100, 100, 1500, 250]); 

for sti_i = 1:size(sort_edge_stdResult_value,2)
    subplot(1,size(sort_edge_stdResult_value,2),sti_i)
    try
        curr_sort_all_neuron_stdResult_value = sort_all_neuron_stdResult_value{sti_i};
        curr_sort_all_neuron_stdResult_value = curr_sort_all_neuron_stdResult_value(curr_sort_all_neuron_stdResult_value > 0);
        if isempty(curr_sort_all_neuron_stdResult_value)
            warning('Stimulus %d: Neuron data is empty or all non-positive after filtering.', sti_i);
            title(sprintf('Stimulus %d (No Neuron Data)', sti_i), 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');
            continue; 
        end
        [f_neuron, x_neuron] = ksdensity(curr_sort_all_neuron_stdResult_value, 'Support', 'positive'); % 使用核密度估计
        if isempty(f_neuron) || isempty(x_neuron) || all(isnan(f_neuron)) || all(isnan(x_neuron))
             warning('Stimulus %d: ksdensity for neuron returned empty or NaN.', sti_i);
             title(sprintf('Stimulus %d (Neuron ksdensity error)', sti_i), 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');
             continue;
        end
        f_neuron = f_neuron / trapz(x_neuron, f_neuron); % 确保密度总和为1

        % 填充区域（在曲线下方着色）
        h_fill_neuron = fill([x_neuron, fliplr(x_neuron)], [f_neuron, zeros(size(f_neuron))], [0.8 0.8 1], 'FaceAlpha', 0.5, 'EdgeColor', 'none');hold on;
        plot(x_neuron, f_neuron, 'DisplayName', 'Neuron Data', 'LineWidth', 0.8, 'LineStyle', '--','Color', color_scheme_npg(3, :));hold on;


        % 计算边密度估计
        current_edge_data = sort_edge_stdResult_value{sti_i};
        if process_nozeros == 1
            aa = nonzeros(current_edge_data);
            if isempty(aa)
                warning('Stimulus %d: Edge data is all zeros or empty after nonzeros().', sti_i);
                title(sprintf('Stimulus %d (No Edge Data)', sti_i), 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');
                set(gca, 'XScale', 'log'); 
                xlim_min = 0.1; 
                if min(x_neuron) < xlim_min && min(x_neuron) > 0
                     xlim_min = max(0.001, min(x_neuron)*0.5); 
                end
                xlim_max_val = 200;
                if max(x_neuron) > xlim_max_val
                    xlim_max_val = max(x_neuron) * 1.1;
                end
                xlim([xlim_min xlim_max_val]); 
                xlabel('Trial Variance (log scale)', 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold'); 
                ylabel('Probability Density', 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold'); 
                set(gca, 'FontSize', 8, 'FontName', 'Arial'); 
                legend([h_fill_neuron], {'All Neuron'}, 'Location', 'best', 'FontSize', 6);
                continue;
            end
            num_to_select = min(select_neuron_num, length(aa));
            if num_to_select == 0 % 如果筛选后还是没有数据
                 warning('Stimulus %d: Edge data (nonzeros) has zero elements to select.', sti_i);
                 title(sprintf('Stimulus %d (No Edge Data to Select)', sti_i), 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');
                set(gca, 'XScale', 'log'); 
                xlim_min = 0.1; 
                if min(x_neuron) < xlim_min && min(x_neuron) > 0
                     xlim_min = max(0.001, min(x_neuron)*0.5);
                end
                xlim_max_val = 200;
                 if max(x_neuron) > xlim_max_val
                    xlim_max_val = max(x_neuron) * 1.1;
                end
                xlim([xlim_min xlim_max_val]); 
                xlabel('Trial Variance (log scale)', 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');
                ylabel('Probability Density', 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');
                set(gca, 'FontSize', 8, 'FontName', 'Arial');
                legend([h_fill_neuron], {'All Neuron'}, 'Location', 'best', 'FontSize', 6);
                continue;
            end
            % aa_selected = aa(1:num_to_select);
            aa_selected = aa;
            [f_edge, x_edge] = ksdensity(aa_selected, 'Support', 'positive');
        else
            current_edge_data = current_edge_data(current_edge_data > 0); % 确保正值
            if isempty(current_edge_data)
                warning('Stimulus %d: Edge data is empty or all non-positive after filtering (process_nozeros=0).', sti_i);
                title(sprintf('Stimulus %d (No Edge Data)', sti_i), 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');
                 set(gca, 'XScale', 'log'); 
                xlim_min = 0.1; 
                if min(x_neuron) < xlim_min && min(x_neuron) > 0
                     xlim_min = max(0.001, min(x_neuron)*0.5);
                end
                xlim_max_val = 200;
                 if max(x_neuron) > xlim_max_val
                    xlim_max_val = max(x_neuron) * 1.1;
                end
                xlim([xlim_min xlim_max_val]); 
                xlabel('Trial Variance (log scale)', 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');
                ylabel('Probability Density', 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');
                set(gca, 'FontSize', 8, 'FontName', 'Arial');
                legend([h_fill_neuron], {'All Neuron'}, 'Location', 'best', 'FontSize', 6);
                continue;
            end
            [f_edge, x_edge] = ksdensity(current_edge_data, 'Support', 'positive'); % 使用核密度估计
        end

        if isempty(f_edge) || isempty(x_edge) || all(isnan(f_edge)) || all(isnan(x_edge))
             warning('Stimulus %d: ksdensity for edge returned empty or NaN.', sti_i);
             title(sprintf('Stimulus %d (Edge ksdensity error)', sti_i), 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');
             % 只绘制神经元数据
             set(gca, 'XScale', 'log');
             xlim_min = 0.1; 
             if min(x_neuron) < xlim_min && min(x_neuron) > 0
                 xlim_min = max(0.001, min(x_neuron)*0.5);
             end
             xlim_max_val = 200;
             if max(x_neuron) > xlim_max_val
                xlim_max_val = max(x_neuron) * 1.1;
             end
             xlim([xlim_min xlim_max_val]); 
             xlabel('Trial Variance (log scale)', 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');
             ylabel('Probability Density', 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');
             set(gca, 'FontSize', 8, 'FontName', 'Arial');
             legend([h_fill_neuron], {'All Neuron'}, 'Location', 'best', 'FontSize', 6);
             continue;
        end
        f_edge = f_edge / trapz(x_edge, f_edge); 

        h_fill_edge = fill([x_edge, fliplr(x_edge)], [f_edge, zeros(size(f_edge))],  [1 0.85 0.6], 'FaceAlpha', 0.5, 'EdgeColor', 'none');hold on;
        plot(x_edge, f_edge, 'DisplayName', 'Edge Data', 'LineWidth', 0.8,  'Color', color_scheme_npg(2, :)); 

        set(gca, 'XScale', 'log'); 

        all_x_values = [x_neuron, x_edge];
        all_x_values_positive = all_x_values(all_x_values > 0);
        if isempty(all_x_values_positive)
            xlim_min_val = 0.01;
            xlim_max_val = 200;   
        else
            xlim_min_val = max(0.001, min(all_x_values_positive) * 0.5); 
            xlim_max_val = max(all_x_values_positive) * 1.5;      
            if xlim_max_val < 10 
                xlim_max_val = 10;
            end
            if xlim_min_val >= xlim_max_val 
                xlim_min_val = xlim_max_val * 0.01;
            end
        end

        xlim([xlim_min_val xlim_max_val]); 

        xlabel('Trial Variance (log scale)', 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');   
        ylabel('Probability Density', 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');    
        set(gca, 'FontSize', 8, 'FontName', 'Arial');            

        legend([h_fill_neuron, h_fill_edge], {'All Neuron', 'All Edge'}, 'Location', 'best', 'FontSize', 6);
        title(sprintf('Stimulus %d', sti_i), 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');
    catch ME 
        warning('Error processing subplot %d: %s', sti_i, ME.message);
        fprintf('Error in file %s, line %d\n', ME.stack(1).file, ME.stack(1).line);
        title(sprintf('Stimulus %d (Error)', sti_i), 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');
    end
end
if ~exist(figure_path, 'dir')
   mkdir(figure_path);
end

savefig([figure_path sprintf('trial_variability_all_neuron&edge_PDF_statistic_logX.fig')]); 
exportgraphics(gcf,[figure_path sprintf('trial_variability_all_neuron&edge_PDF_statistic_logX'),'.png'],'Resolution',300)





figure('Position', [100, 100, 600, 250]); 

for sti_i = 1:size(compare_sort_edge_stdResult_value,2)
    subplot(1,size(compare_sort_edge_stdResult_value,2),sti_i)
    try
        % 计算神经元密度估计
        curr_sort_all_neuron_stdResult_value = compare_sort_all_neuron_stdResult_value{sti_i};
        curr_sort_all_neuron_stdResult_value = curr_sort_all_neuron_stdResult_value(curr_sort_all_neuron_stdResult_value > 0);
        if isempty(curr_sort_all_neuron_stdResult_value)
            warning('Stimulus %d: Neuron data is empty or all non-positive after filtering.', sti_i);
            title(sprintf('Stimulus %d (No Neuron Data)', sti_i), 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');
            continue; 
        end
        [f_neuron, x_neuron] = ksdensity(curr_sort_all_neuron_stdResult_value, 'Support', 'positive'); 
        if isempty(f_neuron) || isempty(x_neuron) || all(isnan(f_neuron)) || all(isnan(x_neuron))
             warning('Stimulus %d: ksdensity for neuron returned empty or NaN.', sti_i);
             title(sprintf('Stimulus %d (Neuron ksdensity error)', sti_i), 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');
             continue;
        end
        f_neuron = f_neuron / trapz(x_neuron, f_neuron); 

        % 填充区域（在曲线下方着色）
        h_fill_neuron = fill([x_neuron, fliplr(x_neuron)], [f_neuron, zeros(size(f_neuron))], [0.8 0.8 1], 'FaceAlpha', 0.5, 'EdgeColor', 'none');hold on;
        plot(x_neuron, f_neuron, 'DisplayName', 'Neuron Data', 'LineWidth', 0.8, 'LineStyle', '--','Color', color_scheme_npg(3, :));hold on;

        % 计算边密度估计
        current_edge_data = compare_sort_edge_stdResult_value{sti_i};
        if process_nozeros == 1
            aa = nonzeros(current_edge_data);
            if isempty(aa)
                warning('Stimulus %d: Edge data is all zeros or empty after nonzeros().', sti_i);
                title(sprintf('Stimulus %d (No Edge Data)', sti_i), 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');
                set(gca, 'XScale', 'log');
                xlim_min = 0.1; 
                if min(x_neuron) < xlim_min && min(x_neuron) > 0
                     xlim_min = max(0.001, min(x_neuron)*0.5); 
                end
                xlim_max_val = 200;
                if max(x_neuron) > xlim_max_val
                    xlim_max_val = max(x_neuron) * 1.1;
                end
                xlim([xlim_min xlim_max_val]); 
                xlabel('Trial Variance (log scale)', 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold'); 
                ylabel('Probability Density', 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold'); 
                set(gca, 'FontSize', 8, 'FontName', 'Arial'); 
                legend([h_fill_neuron], {'All Neuron'}, 'Location', 'best', 'FontSize', 6);
                continue;
            end
            num_to_select = min(select_neuron_num, length(aa));
            if num_to_select == 0 
                 warning('Stimulus %d: Edge data (nonzeros) has zero elements to select.', sti_i);
                 title(sprintf('Stimulus %d (No Edge Data to Select)', sti_i), 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');
                set(gca, 'XScale', 'log'); 
                xlim_min = 0.1; 
                if min(x_neuron) < xlim_min && min(x_neuron) > 0
                     xlim_min = max(0.001, min(x_neuron)*0.5);
                end
                xlim_max_val = 200;
                 if max(x_neuron) > xlim_max_val
                    xlim_max_val = max(x_neuron) * 1.1;
                end
                xlim([xlim_min xlim_max_val]); 
                xlabel('Trial Variance (log scale)', 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');
                ylabel('Probability Density', 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');
                set(gca, 'FontSize', 8, 'FontName', 'Arial');
                legend([h_fill_neuron], {'All Neuron'}, 'Location', 'best', 'FontSize', 6);
                continue;
            end
            % aa_selected = aa(1:num_to_select);
            aa_selected = aa;
            [f_edge, x_edge] = ksdensity(aa_selected, 'Support', 'positive');
        else
            current_edge_data = current_edge_data(current_edge_data > 0); % 确保正值
            if isempty(current_edge_data)
                warning('Stimulus %d: Edge data is empty or all non-positive after filtering (process_nozeros=0).', sti_i);
                title(sprintf('Stimulus %d (No Edge Data)', sti_i), 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');
                 set(gca, 'XScale', 'log'); 
                xlim_min = 0.1; 
                if min(x_neuron) < xlim_min && min(x_neuron) > 0
                     xlim_min = max(0.001, min(x_neuron)*0.5);
                end
                xlim_max_val = 200;
                 if max(x_neuron) > xlim_max_val
                    xlim_max_val = max(x_neuron) * 1.1;
                end
                xlim([xlim_min xlim_max_val]); 
                xlabel('Trial Variance (log scale)', 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');
                ylabel('Probability Density', 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');
                set(gca, 'FontSize', 8, 'FontName', 'Arial');
                legend([h_fill_neuron], {'All Neuron'}, 'Location', 'best', 'FontSize', 6);
                continue;
            end
            [f_edge, x_edge] = ksdensity(current_edge_data, 'Support', 'positive'); % 使用核密度估计
        end

        if isempty(f_edge) || isempty(x_edge) || all(isnan(f_edge)) || all(isnan(x_edge))
             warning('Stimulus %d: ksdensity for edge returned empty or NaN.', sti_i);
             title(sprintf('Stimulus %d (Edge ksdensity error)', sti_i), 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');
             set(gca, 'XScale', 'log');
             xlim_min = 0.1; 
             if min(x_neuron) < xlim_min && min(x_neuron) > 0
                 xlim_min = max(0.001, min(x_neuron)*0.5);
             end
             xlim_max_val = 200;
             if max(x_neuron) > xlim_max_val
                xlim_max_val = max(x_neuron) * 1.1;
             end
             xlim([xlim_min xlim_max_val]); 
             xlabel('Trial Variance (log scale)', 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');
             ylabel('Probability Density', 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');
             set(gca, 'FontSize', 8, 'FontName', 'Arial');
             legend([h_fill_neuron], {'All Neuron'}, 'Location', 'best', 'FontSize', 6);
             continue;
        end
        f_edge = f_edge / trapz(x_edge, f_edge); % 确保密度总和为1

        h_fill_edge = fill([x_edge, fliplr(x_edge)], [f_edge, zeros(size(f_edge))],  [1 0.85 0.6], 'FaceAlpha', 0.5, 'EdgeColor', 'none');hold on;
        plot(x_edge, f_edge, 'DisplayName', 'Edge Data', 'LineWidth', 0.8,  'Color', color_scheme_npg(2, :)); 

        set(gca, 'XScale', 'log'); 

        all_x_values = [x_neuron, x_edge];
        all_x_values_positive = all_x_values(all_x_values > 0);
        if isempty(all_x_values_positive)
            xlim_min_val = 0.01; 
            xlim_max_val = 200;   
        else
            xlim_min_val = max(0.001, min(all_x_values_positive) * 0.5); 
            xlim_max_val = max(all_x_values_positive) * 1.5;       
            if xlim_max_val < 10 
                xlim_max_val = 10;
            end
            if xlim_min_val >= xlim_max_val 
                xlim_min_val = xlim_max_val * 0.01;
            end
        end

        xlim([xlim_min_val xlim_max_val]); 

        xlabel('Trial Variance (log scale)', 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');    
        ylabel('Probability Density', 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');    
        set(gca, 'FontSize', 8, 'FontName', 'Arial');            

        legend([h_fill_neuron, h_fill_edge], {'All Neuron', 'All Edge'}, 'Location', 'best', 'FontSize', 6);
        title(sprintf('Stimulus %d', sti_i), 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');
    catch ME 
        warning('Error processing subplot %d: %s', sti_i, ME.message);
        fprintf('Error in file %s, line %d\n', ME.stack(1).file, ME.stack(1).line);
        title(sprintf('Stimulus %d (Error)', sti_i), 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');
    end
end
if ~exist(figure_path, 'dir')
   mkdir(figure_path);
end

savefig([figure_path sprintf('trial_variability_all_neuron&edge_PDF_statistic_logX_average_sti.fig')]); 
exportgraphics(gcf,[figure_path sprintf('trial_variability_all_neuron&edge_PDF_statistic_logX_average_sti'),'.png'],'Resolution',300)



% % select neuron & all edge top TVV
figure('Position', [100, 100, 1500, 250]); 
bin_min = 0;
bin_max = 80;
num_bins = 50; % 您可以根据需要调整分箱数量
bin_edges = linspace(bin_min, bin_max, num_bins + 1);

for sti_i = 1:size(sort_edge_stdResult_value,2)
    subplot(1,size(sort_edge_stdResult_value,2),sti_i)
    hold on; 

    try
        % --- 神经元数据直方图 ---
        curr_sort_all_neuron_stdResult_value = sort_all_neuron_stdResult_value{sti_i};

        neuron_face_color = color_scheme_npg(3, :); 
        % neuron_face_color = [0.8 0.8 1]; 

        h_hist_neuron = histogram(curr_sort_all_neuron_stdResult_value, ...
                                  'BinEdges', bin_edges, ...
                                  'Normalization', 'count', ... 
                                  'FaceColor', neuron_face_color, ...
                                  'FaceAlpha', 0.6, ... 
                                  'EdgeColor', 'none'); 

        % --- 边缘数据直方图 ---
        edge_data_for_hist = [];
        if process_nozeros == 1
            aa = nonzeros(sort_edge_stdResult_value{sti_i});
            if length(aa) >= select_neuron_num
                 % edge_data_for_hist = aa(1:select_neuron_num);
                 edge_data_for_hist = aa;
            else
                 edge_data_for_hist = aa; 
                 if ~isempty(aa)
                    fprintf('警告: Stimulus %d, 非零边缘数据点 (%d) 少于 select_neuron_num (%d)。使用所有非零点。\n', sti_i, length(aa), select_neuron_num);
                 else
                    fprintf('警告: Stimulus %d, 非零边缘数据点为空。\n', sti_i);
                 end
            end
        else
            edge_data_for_hist = sort_edge_stdResult_value{sti_i};
        end
        
        edge_face_color = color_scheme_npg(2, :); 
        % edge_face_color = [1 0.85 0.6]; 

        if ~isempty(edge_data_for_hist) 
            h_hist_edge = histogram(edge_data_for_hist, ...
                                    'BinEdges', bin_edges, ...
                                    'Normalization', 'count', ... 
                                    'FaceColor', edge_face_color, ...
                                    'FaceAlpha', 0.6, ... 
                                    'EdgeColor', 'none'); 
        else
            h_hist_edge = []; 
        end

        % --- 坐标轴和标签 ---
        xlim([bin_min bin_max]); 
        xlabel('Trial Variance', 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');
        ylabel('Count', 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold'); 
        set(gca, 'FontSize', 8, 'FontName', 'Arial');

        % --- 添加图例与样式 ---
        legend_handles = [];
        legend_labels = {};
        if ~isempty(h_hist_neuron) && h_hist_neuron.NumBins > 0 
            legend_handles(end+1) = h_hist_neuron;
            legend_labels{end+1} = 'All Neuron';
        end
        if ~isempty(h_hist_edge) && h_hist_edge.NumBins > 0 
            legend_handles(end+1) = h_hist_edge;
            legend_labels{end+1} = 'All Edge';
        end
        
        if ~isempty(legend_handles)
            legend(legend_handles, legend_labels, 'Location', 'best', 'FontSize', 6, 'Box', 'off');
        end
        
        title(sprintf('Stimulus %d', sti_i), 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');
        
    catch ME
        fprintf('处理 Stimulus %d 时发生错误: %s\n', sti_i, ME.message);
    end
    hold off; 
end

try
    savefig_path = fullfile(figure_path, 'trial_variability_all_neuron&edge_HIST_statistic.fig');
    export_path = fullfile(figure_path, 'trial_variability_all_neuron&edge_HIST_statistic.png');
    
    savefig(gcf, savefig_path);
    exportgraphics(gcf, export_path, 'Resolution', 300);
    fprintf('图形已保存到: \n%s\n%s\n', savefig_path, export_path);
catch ME_save
    fprintf('保存图形时发生错误: %s\n', ME_save.message);
end



% % select neuron & all edge top TVV
figure('Position', [100, 100, 600, 250]); 
bin_min = 0;
bin_max = 80;
num_bins = 50; 
bin_edges = linspace(bin_min, bin_max, num_bins + 1);

for sti_i = 1:size(compare_sort_edge_stdResult_value,2)
    subplot(1,size(compare_sort_edge_stdResult_value,2),sti_i)
    hold on; 

    try
        curr_sort_all_neuron_stdResult_value = compare_sort_all_neuron_stdResult_value{sti_i};
        neuron_face_color = color_scheme_npg(3, :); 
        % neuron_face_color = [0.8 0.8 1]; 

        h_hist_neuron = histogram(curr_sort_all_neuron_stdResult_value, ...
                                  'BinEdges', bin_edges, ...
                                  'Normalization', 'count', ... 
                                  'FaceColor', neuron_face_color, ...
                                  'FaceAlpha', 0.6, ... 
                                  'EdgeColor', 'none'); 

        edge_data_for_hist = [];
        if process_nozeros == 1
            aa = nonzeros(compare_sort_edge_stdResult_value{sti_i});
            if length(aa) >= select_neuron_num
                 % edge_data_for_hist = aa(1:select_neuron_num);
                 edge_data_for_hist = aa;
            else
                 edge_data_for_hist = aa; 
                 if ~isempty(aa)
                    fprintf('警告: Stimulus %d, 非零边缘数据点 (%d) 少于 select_neuron_num (%d)。使用所有非零点。\n', sti_i, length(aa), select_neuron_num);
                 else
                    fprintf('警告: Stimulus %d, 非零边缘数据点为空。\n', sti_i);
                 end
            end
        else
            edge_data_for_hist = compare_sort_edge_stdResult_value{sti_i};
        end
        
        edge_face_color = color_scheme_npg(2, :); 

        if ~isempty(edge_data_for_hist) 
            h_hist_edge = histogram(edge_data_for_hist, ...
                                    'BinEdges', bin_edges, ...
                                    'Normalization', 'count', ... 
                                    'FaceColor', edge_face_color, ...
                                    'FaceAlpha', 0.6, ... 
                                    'EdgeColor', 'none'); 
        else
            h_hist_edge = []; 
        end

        xlim([bin_min bin_max]); 
        xlabel('Trial Variance', 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');
        ylabel('Count', 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold'); 
        set(gca, 'FontSize', 8, 'FontName', 'Arial');

        legend_handles = [];
        legend_labels = {};
        if ~isempty(h_hist_neuron) && h_hist_neuron.NumBins > 0 
            legend_handles(end+1) = h_hist_neuron;
            legend_labels{end+1} = 'All Neuron';
        end
        if ~isempty(h_hist_edge) && h_hist_edge.NumBins > 0 
            legend_handles(end+1) = h_hist_edge;
            legend_labels{end+1} = 'All Edge';
        end
        
        if ~isempty(legend_handles)
            legend(legend_handles, legend_labels, 'Location', 'best', 'FontSize', 6, 'Box', 'off');
        end
        
        title(sprintf('Stimulus %d', sti_i), 'FontSize', 8, 'FontName', 'Arial','FontWeight', 'bold');
        
    catch ME
    end
    hold off; 
end

try
    savefig_path = fullfile(figure_path, 'trial_variability_all_neuron&edge_HIST_statistic_average_sti.fig');
    export_path = fullfile(figure_path, 'trial_variability_all_neuron&edge_HIST_statistic_average_sti.png'); 
    savefig(gcf, savefig_path);
    exportgraphics(gcf, export_path, 'Resolution', 300);
catch ME_save
end


