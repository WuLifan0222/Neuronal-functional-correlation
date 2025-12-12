clc; clear; close all; set(0,'defaultfigurecolor','w');
%% generate_neuronal_functional_connectivity
% WLF 20230620
addpath(genpath('D:\PPT_gather\Neuronal functional connectivity\Utils'));

all_start_tim= tic;

%% >>>>>>>>> load color
color_scheme_npg = get_color();

%% >>>>>>>>> Load experiment (load data and stimuli)
folder_path = ['D:\PPT_gather\Neuronal functional connectivity\dataset\'];

param.res_nums_thre = 8; % change here select neurons by response nums % adjust here 1
figure_path = [folder_path sprintf('result_for_thre%d\\',param.res_nums_thre)];
result_for_response_to_stis_path = [figure_path '\'];
mkdir(result_for_response_to_stis_path);
input_C_data = [result_for_response_to_stis_path 'neuron_region_lr_with_select_neurons.mat'];
load(input_C_data);
C = C_mat;

% >>>>>>>> load stimuli
% >>>>>>>>> find different stis
load([folder_path 'demo_visual_stimuli_with_label.mat']);
stimuli = load([folder_path 'demo_visual_stimuli_with_label.mat']);

% set processing param for debug
generate_different_win_corr_result = 1;
generate_different_win_corr_result_for_sub_region = 0;
generate_different_win_corr_result_for_sub_region_only_select = 0;
generate_corr_for_all_neurons = 1;
analysis_all_intensity = 1;
analysis_all_corr_connection = 1;
select_corr_edge_by_statistic_way_for_each_sti = 0;
select_corr_edge_by_trial_variability_way = 1;
top_edge = 200;
top_neuron = 20;

% padding data for tail
C = [C,zeros(size(C,1),40)];
select_neurons_C = [select_neurons_C,zeros(size(select_neurons_C,1),40)];
[select_neurons_C_trials_cell,stimuli_sort,select_neurons_C_sort_trace,before, during, after, sti_labels,sti_num, trials_num] = sort_multi_neurons_trials_matrix_by_sti(select_neurons_C, stimuli, param);

frame_n= size(C,2);
param.show_corr_image = 1;

%% >>>>>>>>> atlas
atlas_bg = importdata('D:\PPT_gather\Neuronal functional connectivity\Utils\atlas\cortical_out_line_resize_5.mat');
atlas = importdata('D:\PPT_gather\Neuronal functional connectivity\Utils\atlas\atlas_top_projection.mat');
clean_top_projection = importdata('D:\PPT_gather\Neuronal functional connectivity\Utils\atlas\clean_top_projection.mat');
resize_ratio = 5;

%% >>>>>> ---------------------------------------------------------------------------------------------------------- generate different win corr result
if generate_different_win_corr_result == 1
    tic;
    corr_param.corr_win = 12;
    corr_param.corr_win_before_num = param.show_sti_before*param.fs/corr_param.corr_win;
    corr_param.corr_win_during_num = param.sti_during_time*param.fs/corr_param.corr_win;
    corr_param.corr_win_after_num = param.show_sti_after*param.fs/corr_param.corr_win;
    corr_param.corr_win_num = size(C,2)-corr_param.corr_win;

    try  
        [select_neurons_C_sort_by_region_corr_matrix, select_neurons_C_sort_by_region_edge_cell] = get_corr_matrix_from_neuron_matrix(select_neurons_C,corr_param);
    catch
        select_neurons_C_sort_by_region_corr_matrix =[];
        select_neurons_C_sort_by_region_edge_cell = {};
    end
    elapsedTime = toc;
    tic
    save([result_for_response_to_stis_path 'neuron_region_lr_with_select_neurons.mat'],...
        'select_neurons_C_sort_by_region_edge_cell',...
        'select_neurons_C_sort_by_region_corr_matrix', ...
        'param','stimuli','-append');

    save([result_for_response_to_stis_path 'data_to_tt.mat'],...
        'C','select_neurons_C','select_neurons_id','select_neurons_C_sort_by_region_edge_cell',...
        'select_neurons_C_sort_by_region_corr_matrix', ...
        'param','stimuli','-v7.3');
    clear select_neurons_C_sort_by_region_edge_cell;
    clear select_neurons_C_sort_by_region_corr_matrix;
    elapsedTime = toc;

    tic;
    % save different win corr
    corr_param = save_different_corr_matrix_by_window_size(C,select_neurons_C,2,param,stimuli,result_for_response_to_stis_path);
    corr_param = save_different_corr_matrix_by_window_size(C,select_neurons_C,4,param,stimuli,result_for_response_to_stis_path);
    corr_param = save_different_corr_matrix_by_window_size(C,select_neurons_C,8,param,stimuli,result_for_response_to_stis_path);
    corr_param = save_different_corr_matrix_by_window_size(C,select_neurons_C,12,param,stimuli,result_for_response_to_stis_path);
    corr_param = save_different_corr_matrix_by_window_size(C,select_neurons_C,16,param,stimuli,result_for_response_to_stis_path);
    corr_param = save_different_corr_matrix_by_window_size(C,select_neurons_C,20,param,stimuli,result_for_response_to_stis_path);
    corr_param = save_different_corr_matrix_by_window_size(C,select_neurons_C,24,param,stimuli,result_for_response_to_stis_path);
    corr_param = save_different_corr_matrix_by_window_size(C,select_neurons_C,36,param,stimuli,result_for_response_to_stis_path);
    corr_param = save_different_corr_matrix_by_window_size(C,select_neurons_C,48,param,stimuli,result_for_response_to_stis_path);
    corr_param = save_different_corr_matrix_by_window_size(C,select_neurons_C,52,param,stimuli,result_for_response_to_stis_path);
    corr_param = save_different_corr_matrix_by_window_size(C,select_neurons_C,60,param,stimuli,result_for_response_to_stis_path);
    corr_param = save_different_corr_matrix_by_window_size(C,select_neurons_C,96,param,stimuli,result_for_response_to_stis_path);
    elapsedTime = toc;

end

%% >>>>>> ---------------------------------------------------------------------------------------------------------- generate different win corr result for sub region 
if generate_different_win_corr_result_for_sub_region == 1

    corr_param.corr_win = 12;
    corr_param.corr_win_before_num = param.show_sti_before*param.fs/corr_param.corr_win;
    corr_param.corr_win_during_num = param.sti_during_time*param.fs/corr_param.corr_win;
    corr_param.corr_win_after_num = param.show_sti_after*param.fs/corr_param.corr_win;
    corr_param.corr_win_num = size(C,2)-corr_param.corr_win;

    all_region_left_neurons_C_sort_by_region_corr_matrix = {};
    all_region_left_neurons_C_sort_by_region_edge_cell = {};
    all_region_right_neurons_C_sort_by_region_corr_matrix = {};
    all_region_right_neurons_C_sort_by_region_edge_cell = {};

    all_region_select_left_neurons_C_sort_by_region_corr_matrix = {};
    all_region_select_left_neurons_C_sort_by_region_edge_cell = {};
    all_region_select_right_neurons_C_sort_by_region_corr_matrix = {};
    all_region_select_right_neurons_C_sort_by_region_edge_cell = {};

    for sub_region_i = 1:size(unique_brain_region,2)
        curr_region_left = [];
        curr_left_neurons_C_sort_by_region_corr_matrix =[];
        curr_left_neurons_C_sort_by_region_edge_cell = {};
        curr_region_left = C_array_left{1,sub_region_i};
        try  
            [curr_left_neurons_C_sort_by_region_corr_matrix, curr_left_neurons_C_sort_by_region_edge_cell] = get_corr_matrix_from_neuron_matrix(curr_region_left,corr_param);
        catch
            curr_left_neurons_C_sort_by_region_corr_matrix =[];
            curr_left_neurons_C_sort_by_region_edge_cell = {};
        end
        all_region_left_neurons_C_sort_by_region_corr_matrix{sub_region_i} = curr_left_neurons_C_sort_by_region_corr_matrix;
        all_region_left_neurons_C_sort_by_region_edge_cell{sub_region_i} = curr_left_neurons_C_sort_by_region_edge_cell;

        curr_region_right = [];
        curr_right_neurons_C_sort_by_region_corr_matrix =[];
        curr_right_neurons_C_sort_by_region_edge_cell = {};
        curr_region_right = C_array_right{1,sub_region_i};
        try  
            [curr_right_neurons_C_sort_by_region_corr_matrix, curr_right_neurons_C_sort_by_region_edge_cell] = get_corr_matrix_from_neuron_matrix(curr_region_right,corr_param);
        catch
            curr_right_neurons_C_sort_by_region_corr_matrix =[];
            curr_right_neurons_C_sort_by_region_edge_cell = {};
        end
        all_region_right_neurons_C_sort_by_region_corr_matrix{sub_region_i} = curr_right_neurons_C_sort_by_region_corr_matrix;
        all_region_right_neurons_C_sort_by_region_edge_cell{sub_region_i} = curr_right_neurons_C_sort_by_region_edge_cell;

        curr_region_select_left = [];
        curr_select_left_neurons_C_sort_by_region_corr_matrix =[];
        curr_select_left_neurons_C_sort_by_region_edge_cell = {};
        curr_region_select_left = select_neurons_sub_left_region_C{1,sub_region_i};
        try  
            [curr_select_left_neurons_C_sort_by_region_corr_matrix, curr_select_left_neurons_C_sort_by_region_edge_cell] = get_corr_matrix_from_neuron_matrix(curr_region_select_left,corr_param);
        catch
            curr_select_left_neurons_C_sort_by_region_corr_matrix =[];
            curr_select_left_neurons_C_sort_by_region_edge_cell = {};
        end
        all_region_select_left_neurons_C_sort_by_region_corr_matrix{sub_region_i} = curr_select_left_neurons_C_sort_by_region_corr_matrix;
        all_region_select_left_neurons_C_sort_by_region_edge_cell{sub_region_i} = curr_select_left_neurons_C_sort_by_region_edge_cell;

        curr_region_select_right = []
        curr_select_right_neurons_C_sort_by_region_corr_matrix =[];
        curr_select_right_neurons_C_sort_by_region_edge_cell = {};
        curr_region_select_right = select_neurons_sub_right_region_C{1,sub_region_i};
        try  
            [curr_select_right_neurons_C_sort_by_region_corr_matrix, curr_select_right_neurons_C_sort_by_region_edge_cell] = get_corr_matrix_from_neuron_matrix(curr_region_select_right,corr_param);
        catch
            curr_select_right_neurons_C_sort_by_region_corr_matrix =[];
            curr_select_right_neurons_C_sort_by_region_edge_cell = {};
        end
        all_region_select_right_neurons_C_sort_by_region_corr_matrix{sub_region_i} = curr_select_right_neurons_C_sort_by_region_corr_matrix;
        all_region_select_right_neurons_C_sort_by_region_edge_cell{sub_region_i} = curr_select_right_neurons_C_sort_by_region_edge_cell;
    end
    save([result_for_response_to_stis_path 'sub_region_corr.mat'],...
        'all_region_left_neurons_C_sort_by_region_edge_cell','all_region_right_neurons_C_sort_by_region_edge_cell',...
        'all_region_left_neurons_C_sort_by_region_corr_matrix','all_region_right_neurons_C_sort_by_region_corr_matrix', ...
        'all_region_select_left_neurons_C_sort_by_region_edge_cell','all_region_select_right_neurons_C_sort_by_region_edge_cell',...
        'all_region_select_left_neurons_C_sort_by_region_corr_matrix','all_region_select_right_neurons_C_sort_by_region_corr_matrix', ...
        'param','stimuli','-v7.3');
end

%% >>>>>> ---------------------------------------------------------------------------------------------------------- generate different win corr result for sub region for select
if generate_different_win_corr_result_for_sub_region_only_select == 1

    corr_param.corr_win = 12;
    corr_param.corr_win_before_num = param.show_sti_before*param.fs/corr_param.corr_win;
    corr_param.corr_win_during_num = param.sti_during_time*param.fs/corr_param.corr_win;
    corr_param.corr_win_after_num = param.show_sti_after*param.fs/corr_param.corr_win;
    corr_param.corr_win_num = size(C,2)-corr_param.corr_win;

    all_region_select_left_neurons_C_sort_by_region_corr_matrix = {};
    all_region_select_left_neurons_C_sort_by_region_edge_cell = {};
    all_region_select_right_neurons_C_sort_by_region_corr_matrix = {};
    all_region_select_right_neurons_C_sort_by_region_edge_cell = {};

    for sub_region_i = 1:size(unique_brain_region,2)
        curr_region_select_left = [];
        curr_select_left_neurons_C_sort_by_region_corr_matrix =[];
        curr_select_left_neurons_C_sort_by_region_edge_cell = {};
        curr_region_select_left = select_neurons_sub_left_region_C{1,sub_region_i};
        try  
            [curr_select_left_neurons_C_sort_by_region_corr_matrix, curr_select_left_neurons_C_sort_by_region_edge_cell] = get_corr_matrix_from_neuron_matrix(curr_region_select_left,corr_param);
        catch
            curr_select_left_neurons_C_sort_by_region_corr_matrix =[];
            curr_select_left_neurons_C_sort_by_region_edge_cell = {};
        end
        all_region_select_left_neurons_C_sort_by_region_corr_matrix{sub_region_i} = curr_select_left_neurons_C_sort_by_region_corr_matrix;
        all_region_select_left_neurons_C_sort_by_region_edge_cell{sub_region_i} = curr_select_left_neurons_C_sort_by_region_edge_cell;

        curr_region_select_right = []
        curr_select_right_neurons_C_sort_by_region_corr_matrix =[];
        curr_select_right_neurons_C_sort_by_region_edge_cell = {};
        curr_region_select_right = select_neurons_sub_right_region_C{1,sub_region_i};
        try  
            [curr_select_right_neurons_C_sort_by_region_corr_matrix, curr_select_right_neurons_C_sort_by_region_edge_cell] = get_corr_matrix_from_neuron_matrix(curr_region_select_right,corr_param);
        catch
            curr_select_right_neurons_C_sort_by_region_corr_matrix =[];
            curr_select_right_neurons_C_sort_by_region_edge_cell = {};
        end
        all_region_select_right_neurons_C_sort_by_region_corr_matrix{sub_region_i} = curr_select_right_neurons_C_sort_by_region_corr_matrix;
        all_region_select_right_neurons_C_sort_by_region_edge_cell{sub_region_i} = curr_select_right_neurons_C_sort_by_region_edge_cell;
    end

    save([result_for_response_to_stis_path 'sub_region_corr_select.mat'],...
        'all_region_select_left_neurons_C_sort_by_region_edge_cell','all_region_select_right_neurons_C_sort_by_region_edge_cell',...
        'all_region_select_left_neurons_C_sort_by_region_corr_matrix','all_region_select_right_neurons_C_sort_by_region_corr_matrix', ...
        'param','stimuli','-v7.3');
end

if generate_corr_for_all_neurons == 1
    tic;
    corr_param.corr_win = 12;
    corr_param.corr_win_before_num = param.show_sti_before*param.fs/corr_param.corr_win;
    corr_param.corr_win_during_num = param.sti_during_time*param.fs/corr_param.corr_win;
    corr_param.corr_win_after_num = param.show_sti_after*param.fs/corr_param.corr_win;
    corr_param.corr_win_num = size(C,2)-corr_param.corr_win;
    try  
        [all_neurons_C_sort_by_region_corr_matrix, all_neurons_C_sort_by_region_edge_cell] = get_corr_matrix_from_neuron_matrix(C,corr_param);
    catch
        all_neurons_C_sort_by_region_corr_matrix =[];
        all_neurons_C_sort_by_region_edge_cell = {};
    end
    elapsedTime = toc;

    tic
    save([result_for_response_to_stis_path 'all_neuron_corr.mat'],...
        'all_neurons_C_sort_by_region_edge_cell',...
        'all_neurons_C_sort_by_region_corr_matrix', ...
        'param','stimuli','-v7.3');
    clear all_neurons_C_sort_by_region_edge_cell;
    clear all_neurons_C_sort_by_region_corr_matrix;
    elapsedTime = toc;
end

% reload data because cleared corr
input_C_data = [result_for_response_to_stis_path 'neuron_region_lr_with_select_neurons_win12.mat'];
load(input_C_data);
C = C_mat;
load([folder_path 'demo_visual_stimuli_with_label.mat']);
stimuli = load([folder_path 'demo_visual_stimuli_with_label.mat']);

stimuli_for_pac = stimuli;
for i = 1 : length(start_edge)
    stimuli_for_pac.stimuli_array_with_label(stimuli.start_edge(i)-param.fs*(param.show_sti_before) : stimuli.end_edge(i)+param.fs*(param.show_sti_after)) = ...
    stimuli.stimuili_label_ind(i);
end

%% >>>>>> ---------------------------------------------------------------------------------------------------------- analysis all intensity
if analysis_all_intensity == 1
    close all

    [coeff_act, score_act,latent_act,tsquared_act,explained_act,mu_single] =pca(select_neurons_C');
    select_neurons_C_pca = score_act(:, 1:3)';
    dists = pdist(select_neurons_C);  
    avgDist = mean(dists);  

    x = select_neurons_C_pca(1,:);
    y = select_neurons_C_pca(2,:);
    z = select_neurons_C_pca(3,:);
    t = 1:min(size(select_neurons_C_pca,2), size(stimuli_for_pac.stimuli_array_with_label,2));

    colors = stimuli_for_pac.stimuli_array_with_label;
    figure, hold on;
    for ii = 1:length(t)-1
        x_patch = [x(ii), x(ii+1), x(ii+1), x(ii)];
        y_patch = [y(ii), y(ii+1), y(ii+1), y(ii)];
        z_patch = [z(ii), z(ii+1), z(ii+1), z(ii)];
        c_patch = [colors(ii), colors(ii+1), colors(ii+1), colors(ii)];  
        patch(x_patch, y_patch, z_patch, c_patch, 'EdgeColor', 'interp', 'LineWidth', 2);
    end
    hold off; 
    colorbar; 
    xlabel('pca1'); 
    ylabel('pca2'); 
    zlabel('pca3'); 
    title('intensity pca change', 'FontSize', 14, 'FontWeight', 'bold'); 


    str_avgDist = sprintf('Average Intra-Class Distance: %.4f', avgDist);
    dim = [.1 .8 .3 .1]; 
    annotation('textbox', dim, 'String', str_avgDist, ...
                 'FitBoxToText', 'on', ...
                 'LineStyle', 'none', ... 
                 'FontSize', 10, ...
                 'FontName', 'Arial'); 
    view(3); 
    grid off; 
    axis off
    savefig([result_for_response_to_stis_path sprintf('intensity pca change')]);
    exportgraphics(gcf,[result_for_response_to_stis_path sprintf('intensity pca change'),'.png'],'Resolution',300)
end

%% >>>>>> ---------------------------------------------------------------------------------------------------------- analysis all corr connection
if analysis_all_corr_connection ==1
    % show all corr pca change
    close all
    k = 1;  
    [idx, Center] = kmeans(select_neurons_C_sort_by_region_corr_matrix, k);  
    [coeff_act, score_act,latent_act,tsquared_act,explained_act,mu_single] =pca(select_neurons_C_sort_by_region_corr_matrix');
    select_neurons_C_sort_by_region_corr_matrix_pca = score_act(:, 1:3)';
    dists = pdist(select_neurons_C_sort_by_region_corr_matrix); 
    avgDist = mean(dists); 

    x = select_neurons_C_sort_by_region_corr_matrix_pca(1,:);
    y = select_neurons_C_sort_by_region_corr_matrix_pca(2,:);
    z = select_neurons_C_sort_by_region_corr_matrix_pca(3,:);
    t = 1:min(size(select_neurons_C_sort_by_region_corr_matrix_pca,2), size(stimuli.stimuli_array_with_label,2));
    colors = stimuli_for_pac.stimuli_array_with_label;
    figure, hold on;
    for ii = 1:length(t)-1
        x_patch = [x(ii), x(ii+1), x(ii+1), x(ii)];
        y_patch = [y(ii), y(ii+1), y(ii+1), y(ii)];
        z_patch = [z(ii), z(ii+1), z(ii+1), z(ii)];
        c_patch = [colors(ii), colors(ii+1), colors(ii+1), colors(ii)]; 
        patch(x_patch, y_patch, z_patch, c_patch, 'EdgeColor', 'interp', 'LineWidth', 2);
    end
    hold off; 
    colorbar; 
    xlabel('pca1');
    ylabel('pca2'); 
    zlabel('pca3'); 
    title('corr pca change', 'FontSize', 14, 'FontWeight', 'bold'); 

    str_avgDist = sprintf('Average Intra-Class Distance: %.4f', avgDist);
    dim = [.1 .8 .3 .1];
    annotation('textbox', dim, 'String', str_avgDist, ...
                 'FitBoxToText', 'on', ...
                 'LineStyle', 'none', ... 
                 'FontSize', 10, ...
                 'FontName', 'Arial'); %,'FontWeight', 'bold'
    view(3); 
    grid off; 
    axis off
    savefig([result_for_response_to_stis_path sprintf('corr pca change')]);
    exportgraphics(gcf,[result_for_response_to_stis_path sprintf('corr pca change'),'.png'],'Resolution',300)
end


%% >>>>>> ---------------------------------------------------------------------------------------------------------- select corr edge by statistic way for each sti
if select_corr_edge_by_statistic_way_for_each_sti == 1
    param.similarity_thre = 400;
    % select stable edge (smilar)
    [edge_selected_by_similarity_thre_ind,edge_selected_by_similarity_thre,edge_select_by_trials_sort_ind,edge_select_by_trials_sort,edge_trials,edge_trial_similar,select_edge_response_by_trial_for_each_sti_condition,...
        each_sti_edge_trial_similar,each_sti_edge_select_by_trials_sort_ind,each_sti_edge_select_by_trials_sort] = ...
        select_edge_by_trials(select_neurons_C_sort_by_region_corr_matrix,param,stimuli);

    % sort for better visual
    [~,~,select_neurons_C_sort_by_region_corr_matrix_visualsort,~, ~, ~, ~,~, ~] = ...
        sort_multi_neurons_trials_matrix_by_sti(select_neurons_C_sort_by_region_corr_matrix, stimuli, param);

    % get corr trials data for show
    [corr_trials_cell,before,during,after,sti_num] = get_corr_trials_matrix_by_sti(select_neurons_C_sort_by_region_edge_cell,stimuli,param);
    C_matNorm = (C_mat - min(C_mat(:))) / (max(C_mat(:)) - min(C_mat(:)));
    [C_matNorm_cell,before,during,after,sti_num] = get_trials_matrix_multi_trace_by_sti((C_matNorm+0.0000001),stimuli,param);
    select_neurons_CNorm = (select_neurons_C - min(select_neurons_C(:))) / (max(select_neurons_C(:)) - min(select_neurons_C(:)));
    [select_neurons_CNorm_cell,before,during,after,sti_num] = get_trials_matrix_multi_trace_by_sti((select_neurons_CNorm+0.0000001),stimuli,param);

    % show result
    for sti_ii = 1:param.sti_num
        curr_select_edge_similarity = each_sti_edge_trial_similar(:,sti_ii);
        % curr_select_edge_id = each_sti_edge_select_by_trials_sort_ind{sti_ii}(1:param.top_edge); % top
        curr_select_edge_id = edge_selected_by_similarity_thre_ind{sti_ii}; % similarity thre
        curr_select_edge_visualsort = select_neurons_C_sort_by_region_corr_matrix_visualsort(curr_select_edge_id,:);
        % >>>>>> show all select edge by deconve way
        deconve_trace_with_label(curr_select_edge_visualsort, stimuli_sort, color_scheme_npg)
        savefig([result_for_response_to_stis_path sprintf('deconve_sti%d_sort_corr_trace_change_constituent',sti_ii)]);
        exportgraphics(gcf,[result_for_response_to_stis_path sprintf('deconve_sti%d_sort_corr_trace_change_constituent',sti_ii),'.png'],'Resolution',300)
        point_num = size(select_neurons_C,1);
        [curr_select_sti_musk] = select_edge_transfer_to_musk(curr_select_edge_id,select_neurons_C_sort_by_region_corr_matrix_visualsort,point_num);
        select_sti_musk{sti_ii} = curr_select_sti_musk;
        figure; imagesc(curr_select_sti_musk)
        savefig([result_for_response_to_stis_path sprintf('select_sti%d_musk',sti_ii)]);
        exportgraphics(gcf,[result_for_response_to_stis_path sprintf('select_sti%d_musk',sti_ii),'.png'],'Resolution',300)
        close all
    end

    save([result_for_response_to_stis_path 'select_edges.mat'], ...
        'edge_select_by_trials_sort','edge_trials','each_sti_edge_trial_similar',...
        'each_sti_edge_select_by_trials_sort_ind','each_sti_edge_select_by_trials_sort','select_sti_musk','edge_selected_by_similarity_thre_ind','edge_selected_by_similarity_thre','-v7.3');

end

if select_corr_edge_by_trial_variability_way == 1
    %% select top stable neurons
    % >>>>>> analysis intensity
    trial_frames = (param.sti_during_time + param.show_sti_after) * param.fs;
    top_N_neurons_C_trace_cell = {};
    top_N_neurons_C_index_cell = {};
    for index = 1:length(select_neurons_C_trials_cell)
        nNeuron = size(select_neurons_C_trials_cell{index},2);
        ori_select_neurons_C_R = reshape(select_neurons_C_trials_cell{index}', nNeuron, trial_frames, []);
        stdResult = zeros(nNeuron,1); 
        for i = 1:nNeuron
            % stdResult(i) = calResNormStd(ori_select_neurons_C_R(i,:, 1:8));
            stdResult(i) = calResNormStd(ori_select_neurons_C_R(i,:, :));
        end
        N = top_neuron;
        [~, sortedIndices] = sort(stdResult);
        minNIndices = sortedIndices(1:N);
        top_N_neurons_C_trace = ori_select_neurons_C_R(minNIndices, :, :);
        top_N_neurons_C_trace_cell{index} = top_N_neurons_C_trace;
        top_N_neurons_C_index_cell{index} = minNIndices;
    end    
    close all 
    save(sprintf('%s/top_%d_neurons_trial_variance.mat', result_for_response_to_stis_path,top_edge), 'top_N_neurons_C_trace_cell','top_N_neurons_C_index_cell','select_neurons_C_sort_trace','-v7.3');
        
    %% select top stable edges
    [select_neurons_C_sort_by_region_corr_matrix_cell,stimuli_sort,select_neurons_C_sort_by_region_corr_matrix_visualsort,~, ~, ~, ~,~, ~] = sort_multi_neurons_trials_matrix_by_sti(select_neurons_C_sort_by_region_corr_matrix, stimuli, param);    
    select_sti_musk = {};
    top_N_matrixs_C_trace_cell = {};
    top_N_matrixs_C_index_cell = {};
    each_sti_edge_select_by_trials_sort = {};
    for index = 1:length(select_neurons_C_sort_by_region_corr_matrix_cell)
        nNeuron = size(select_neurons_C_sort_by_region_corr_matrix_cell{index},2);
        ori_select_matrixs_C_R = reshape(select_neurons_C_sort_by_region_corr_matrix_cell{index}', nNeuron, trial_frames, []);
        stdResult = zeros(nNeuron,1); 
        for i = 1:nNeuron
            % stdResult(i) = calResNormStd(ori_select_matrixs_C_R(i,:, 1:8));
            stdResult(i) = calResNormStd(ori_select_matrixs_C_R(i,:, :));
        end
        % top 100
        N = top_edge;
        [~, sortedIndices] = sort(stdResult);
        minNIndices = sortedIndices(1:N);
        top_N_matrixs_C_trace = ori_select_matrixs_C_R(minNIndices, :, :);
        top_N_C_trace = select_neurons_C_sort_by_region_corr_matrix_cell{index};
        top_N_matrixs_C_trace_cell{index} = top_N_matrixs_C_trace;
        top_N_matrixs_C_index_cell{index} = minNIndices;
        each_sti_edge_select_by_trials_sort{index} = select_neurons_C_sort_by_region_corr_matrix(minNIndices,:);
   
        curr_select_edge_id = minNIndices; % similarity thre
        curr_select_edge_visualsort = select_neurons_C_sort_by_region_corr_matrix_visualsort(curr_select_edge_id,:);

        point_num = size(select_neurons_C,1);
        [curr_select_sti_musk] = select_edge_transfer_to_musk(curr_select_edge_id,select_neurons_C_sort_by_region_corr_matrix_visualsort,point_num);
        select_sti_musk{index} = curr_select_sti_musk;
        figure; imagesc(curr_select_sti_musk)
        savefig([result_for_response_to_stis_path sprintf('select_sti%d_musk',index)]);
        exportgraphics(gcf,[result_for_response_to_stis_path sprintf('select_sti%d_musk',index),'.png'],'Resolution',300)
        close all
    end
    save(sprintf('%s/top_%d_edges_trial_variance.mat', result_for_response_to_stis_path,N),'select_sti_musk','each_sti_edge_select_by_trials_sort','top_N_matrixs_C_trace_cell','top_N_matrixs_C_index_cell','select_neurons_C_sort_by_region_corr_matrix_visualsort','-v7.3');
    
end
