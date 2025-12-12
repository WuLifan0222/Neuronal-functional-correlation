clc; clear; close all; set(0,'defaultfigurecolor','w');
% fig1 adjust
% WLF 20240606

addpath(genpath('D:\PPT_gather\Neuronal functional connectivity\Utils'));

%% ----------------------------------------------------------------------------------------------------------------------------------------------- pre load
% >>>>>>>>> load color
% colormap {light red; light blue; green; blue; light pink; light purple; light green; red; golden; light brown}
color_scheme_npg = get_color();

%% ----------------------------------------------------------------------------------------------------------------------------------------------- show typical trace
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% load data by indicating experiment index
folder_path = ['D:\PPT_gather\Neuronal functional connectivity\dataset\'];
param.res_nums_thre = 8; % change here select neurons by response nums % adjust here 1
figure_path = [folder_path sprintf('result_for_thre%d\\',param.res_nums_thre)];
result_for_response_to_stis_path = [figure_path '\'];
mkdir(result_for_response_to_stis_path);

input_data = [result_for_response_to_stis_path 'index_statistic_all_trial_select_neurons0.95.mat']; 
load(input_data);
try
average_path_length_sequence = select_neurons_average_path_length_sequence;
betweenness_sequence = select_neurons_betweenness_sequence;
closeness_sequence = select_neurons_closeness_sequence;
clustering_sequence = select_neurons_clustering_sequence;
degree_sequence = select_neurons_degree_sequence;
catch
end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% load stimuli
load([folder_path 'demo_visual_stimuli_with_label.mat']);
stimuli = load([folder_path 'demo_visual_stimuli_with_label.mat']);
stimuli_for_edge = load([folder_path 'demo_visual_stimuli_with_label.mat']);
stimuli_for_edge.start_edge = stimuli_for_edge.start_edge - 12;
stimuli_for_edge.end_edge = stimuli_for_edge.end_edge;

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

param.fs = 4;
param.sti_during_time = 10; %change here according different sti
param.show_sti_before = 3; %change here according different sti
param.show_sti_after = 7; %change here according different sti
param.delay_after = 0;
param.before = param.show_sti_before * param.fs;
param.during = param.sti_during_time * param.fs;
param.after = param.show_sti_after * param.fs;

param_for_edge = param;
param_for_edge.sti_during_time = 13; %change here according different sti
param_for_edge.show_sti_before = 0.1; %change here according different sti
param_for_edge.show_sti_after = 7; %change here according different sti


%% show  result
curr_average_path_length_sequence = [];
curr_betweenness_sequence = [];
curr_closeness_sequence = [];
curr_clustering_sequence = [];
curr_degree_sequence = [];
n = 0;

% 计算子图的行数和列数
num_sti = size(average_path_length_sequence, 1);
num_sessions = size(average_path_length_sequence, 2);
nrows = num_sti; % 每个 Sti 一行
ncols = num_sessions; % 每个 Session 一列


figure('Color', 'w','Position', [100, 100, 200, 800]); 
for sti_i = 1:num_sti
    for session_i = 1:num_sessions
        subplot(nrows, ncols, (sti_i - 1) * ncols + session_i); 

        hold on; 

        for trial_i = 1:size(average_path_length_sequence, 3)
            for frame_i = 1:size(average_path_length_sequence, 4)
                curr_average_path_length_sequence(frame_i) = average_path_length_sequence(sti_i, session_i, trial_i, frame_i);
                curr_clustering_sequence(frame_i) = clustering_sequence(sti_i, session_i, trial_i, frame_i);
            end

            xx = 1:size(curr_average_path_length_sequence,2);
            scatter(xx, curr_average_path_length_sequence, ...
            20, color_scheme_npg(9,:), 'filled', ... 
            'LineWidth', 1.5, ... 
            'Marker', 'o');
            plot(curr_clustering_sequence, 'LineWidth', 1.5, 'Color', color_scheme_npg(3,:));

            plot([param.fs * param_for_edge.show_sti_before, param.fs * param_for_edge.show_sti_before], [0, 3], '--r', 'LineWidth', 0.5);
            plot([param.fs * (param_for_edge.show_sti_before + param_for_edge.sti_during_time), param.fs * (param_for_edge.show_sti_before + param_for_edge.sti_during_time)], [0, 3], '--r', 'LineWidth', 0.5);
        end

        xticks = get(gca, 'XTick'); 
        xticklabels = xticks / param.fs - param.show_sti_before; 
        set(gca, 'XTick', xticks, 'XTickLabel', xticklabels); 

        title(sprintf('Sti %d, Session %d', sti_i, session_i), 'FontSize', 12);
        xlabel('Time (s)', 'FontSize', 10); % 修改x轴标签为‘Time (s)’
        ylabel('Values', 'FontSize', 10);
        grid off;
        set(gca, 'FontSize', 10);

        hold off; 
    end
end

subplot(nrows, ncols, nrows*ncols); 
set(gca, 'LooseInset', get(gca, 'TightInset'));
savefig([result_for_response_to_stis_path sprintf('average_path_length_sequence')]);
exportgraphics(gcf, [result_for_response_to_stis_path sprintf('average_path_length_sequence.png')], 'Resolution', 300);



% % >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> fig in sub plot

for sti_i = 1:num_sti
    figure('Color', 'w','Position', [100, 100, 400, 300]); 
    for session_i = 1:num_sessions
        subplot(2, 1, 1); 
        hold on; 
        for trial_i = 1:size(average_path_length_sequence, 3)
            for frame_i = 1:size(average_path_length_sequence, 4)
                curr_average_path_length_sequence(frame_i) = average_path_length_sequence(sti_i, session_i, trial_i, frame_i);
                curr_clustering_sequence(frame_i) = clustering_sequence(sti_i, session_i, trial_i, frame_i);
            end

            xx = 1:size(curr_average_path_length_sequence,2);
            scatter(xx, curr_average_path_length_sequence, ...
            20, color_scheme_npg(9,:), 'filled', ... 
            'LineWidth', 1.5, ... 
            'Marker', 'o'); hold on;

            plot([param.fs * param_for_edge.show_sti_before, param.fs * param_for_edge.show_sti_before], [0, 6], '--r', 'LineWidth', 0.5);hold on;
            plot([param.fs * (param_for_edge.show_sti_before + param_for_edge.sti_during_time), param.fs * (param_for_edge.show_sti_before + param_for_edge.sti_during_time)], [0, 6], '--r', 'LineWidth', 0.5);

        end
        xticks = get(gca, 'XTick'); 
        xticklabels = xticks / param.fs - param.show_sti_before; 
        set(gca, 'XTick', xticks, 'XTickLabel', xticklabels); 

        xlim([0,80]);ylim([0,6])
        title(sprintf('Sti %d, Session %d', sti_i, session_i), 'FontSize', 12);
        xlabel('Time (s)', 'FontSize', 10); % 修改x轴标签为‘Time (s)’
        ylabel('Average Path Length', 'FontSize', 10);
        % grid off;
        % axis off
        box on
        set(gca, 'FontSize', 10);

        subplot(2, 1, 2); 
        for trial_i = 1:size(average_path_length_sequence, 3)
            plot(curr_clustering_sequence, 'LineWidth', 1.5, 'Color', color_scheme_npg(3,:)); hold on;
            plot([param.fs * param_for_edge.show_sti_before, param.fs * param_for_edge.show_sti_before], [0, 1], '--r', 'LineWidth', 0.5);hold on;
            plot([param.fs * (param_for_edge.show_sti_before + param_for_edge.sti_during_time), param.fs * (param_for_edge.show_sti_before + param_for_edge.sti_during_time)], [0, 1], '--r', 'LineWidth', 0.5);
        end 


        xticks = get(gca, 'XTick'); 
        xticklabels = xticks / param.fs - param.show_sti_before; 
        set(gca, 'XTick', xticks, 'XTickLabel', xticklabels); 
        xlim([0,80]);ylim([0,1])
        title(sprintf('Sti %d, Session %d', sti_i, session_i), 'FontSize', 12);
        xlabel('Time (s)', 'FontSize', 10); % 修改x轴标签为‘Time (s)’
        ylabel('Clustering Coefficient', 'FontSize', 10);
        % grid off;
        % axis off
        box on;
        set(gca, 'FontSize', 10);

        hold off; 
    end
    set(gca, 'LooseInset', get(gca, 'TightInset'));
    savefig([result_for_response_to_stis_path sprintf('sti%d_average_path_length_sequence',sti_i)]);
    exportgraphics(gcf, [result_for_response_to_stis_path sprintf('sti%d_average_path_length_sequence.png',sti_i)], 'Resolution', 300);

end



figure('Color', 'w','Position', [100, 100, 200, 800]); 

curr_select_neurons_all_degree_sum_all_trial = {};
curr_select_neurons_all_degree_sum_all_trial_mean = {};
for sti_i = 1:num_sti
% for sti_i = 2
    for session_i = 1:num_sessions
        subplot(nrows, ncols, (sti_i - 1) * ncols + session_i); 
        hold on; 
        for trial_i = 1:size(select_neurons_all_degree, 3)
            for frame_i = 1:size(select_neurons_all_degree, 4)
                curr_select_neurons_all_degree(:,frame_i) = select_neurons_all_degree{sti_i, session_i, trial_i, frame_i};
                curr_select_neurons_all_clustering(frame_i,:) = select_neurons_all_clustering{sti_i, session_i, trial_i, frame_i};
            end

            plot(sum(curr_select_neurons_all_degree,1), 'LineWidth', 1.5, 'Color', color_scheme_npg(9,:)); 

            plot(curr_select_neurons_all_clustering, 'LineWidth', 1.5, 'Color', color_scheme_npg(3,:));
            aa = max(max(sum(curr_select_neurons_all_degree,1)));

            plot([param.fs * param_for_edge.show_sti_before, param.fs * param_for_edge.show_sti_before], [0, aa], '--r', 'LineWidth', 0.5);
            plot([param.fs * (param_for_edge.show_sti_before + param_for_edge.sti_during_time), param.fs * (param_for_edge.show_sti_before + param_for_edge.sti_during_time)], [0, aa], '--r', 'LineWidth', 0.5);
            curr_select_neurons_all_degree_sum_all_trial0(trial_i,:) = sum(curr_select_neurons_all_degree,1);
            
        end
        curr_select_neurons_all_degree_sum_all_trial{sti_i} = curr_select_neurons_all_degree_sum_all_trial0;
        curr_select_neurons_all_degree_sum_all_trial_mean{sti_i} = mean(curr_select_neurons_all_degree_sum_all_trial0,1);

        xticks = get(gca, 'XTick'); 
        xticklabels = xticks / param.fs - param.show_sti_before;
        set(gca, 'XTick', xticks, 'XTickLabel', xticklabels); 

        title(sprintf('Sti %d, Session %d', sti_i, session_i), 'FontSize', 12);
        xlabel('Time (s)', 'FontSize', 10); 
        ylabel('Values', 'FontSize', 10);
        grid off; 
        set(gca, 'FontSize', 10); 

        hold off; 
    end
end
savefig([result_for_response_to_stis_path sprintf('degree_sequence')]);
exportgraphics(gcf, [result_for_response_to_stis_path sprintf('degree_sequence.png')], 'Resolution', 300);



figure('Color', 'w','Position', [100, 100, 200, 800]); 

for sti_i = 1:num_sti
    for session_i = 1:num_sessions
        subplot(nrows, ncols, (sti_i - 1) * ncols + session_i); 
        hold on; 
        for trial_i = 1:size(select_neurons_all_degree, 3)
        % for trial_i = 1
            for frame_i = 1:size(select_neurons_all_degree, 4)
                % 从数据提取值
                curr_select_neurons_all_degree(:,frame_i) = select_neurons_all_degree{sti_i, session_i, trial_i, frame_i};
                curr_select_neurons_all_clustering(frame_i,:) = select_neurons_all_clustering{sti_i, session_i, trial_i, frame_i};
            end

            % 绘制各个序列的曲线
            plot(sum(curr_select_neurons_all_degree,1), 'LineWidth', 1.5, 'Color', color_scheme_npg(9,:)); 

            plot(curr_select_neurons_all_clustering, 'LineWidth', 1.5, 'Color', color_scheme_npg(3,:));
            aa = max(max(sum(curr_select_neurons_all_degree,1)));
            plot([param.fs * param_for_edge.show_sti_before, param.fs * param_for_edge.show_sti_before], [0, aa], '--r', 'LineWidth', 0.5);
            plot([param.fs * (param_for_edge.show_sti_before + param_for_edge.sti_during_time), param.fs * (param_for_edge.show_sti_before + param_for_edge.sti_during_time)], [0, aa], '--r', 'LineWidth', 0.5);
        
        end
        plot(curr_select_neurons_all_degree_sum_all_trial_mean{sti_i}, 'LineWidth', 3, 'Color', color_scheme_npg(1,:)); 

        % 调整x轴以显示秒数
        xticks = get(gca, 'XTick'); 
        xticklabels = xticks / param.fs - param.show_sti_before; 
        set(gca, 'XTick', xticks, 'XTickLabel', xticklabels); 

        title(sprintf('Sti %d, Session %d', sti_i, session_i), 'FontSize', 12);
        xlabel('Time (s)', 'FontSize', 10); % 修改x轴标签为‘Time (s)’
        ylabel('Values', 'FontSize', 10);
        grid off; 
        set(gca, 'FontSize', 10); 

        hold off; 
    end
end
savefig([result_for_response_to_stis_path sprintf('degree_sequence')]);
exportgraphics(gcf, [result_for_response_to_stis_path sprintf('degree_sequence.png')], 'Resolution', 300);


figure('Color', 'w','Position', [100, 100, 200, 800]); 
for sti_i = 1:num_sti
    for session_i = 1:num_sessions
        subplot(nrows, ncols, (sti_i - 1) * ncols + session_i); 
        hold on; 

        plot_boundedline_for_trials_matrix(curr_select_neurons_all_degree_sum_all_trial{sti_i}, color_scheme_npg(1,:)); 
        aa = max(max(mean(curr_select_neurons_all_degree,1),1));
        plot([param.fs * param_for_edge.show_sti_before, param.fs * param_for_edge.show_sti_before], [0, aa], '--r', 'LineWidth', 0.5);
        plot([param.fs * (param_for_edge.show_sti_before + param_for_edge.sti_during_time), param.fs * (param_for_edge.show_sti_before + param_for_edge.sti_during_time)], [0, aa], '--r', 'LineWidth', 0.5);
        
        % 调整x轴以显示秒数
        xticks = get(gca, 'XTick'); 
        xticklabels = xticks / param.fs - param.show_sti_before; 
        set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);

        title(sprintf('Sti %d, Session %d', sti_i, session_i), 'FontSize', 12);
        xlabel('Time (s)', 'FontSize', 10); 
        ylabel('degree', 'FontSize', 10);
        grid off; 
        set(gca, 'FontSize', 10); 

        hold off;
    end
end
savefig([result_for_response_to_stis_path sprintf('degree_sequence_boundline')]);
exportgraphics(gcf, [result_for_response_to_stis_path sprintf('degree_sequence_boundline.png')], 'Resolution', 300);

save([result_for_response_to_stis_path sprintf('select_neurons_all_degree_sum_all_trial.mat')],...
    'curr_select_neurons_all_degree_sum_all_trial_mean',...
    'curr_select_neurons_all_degree_sum_all_trial', ...
    'param','stimuli','-v7.3');



figure('Color', 'w','Position', [100, 100, 150, 800]); 
for sti_i = 1:num_sti
    for session_i = 1:num_sessions
        subplot(nrows, ncols, (sti_i - 1) * ncols + session_i); 
        hold on; 
        for trial_i = 1
            for frame_i = 1:size(select_neurons_all_degree, 4)
                curr_select_neurons_all_degree(:,frame_i) = select_neurons_all_degree{sti_i, session_i, trial_i, frame_i};
                curr_select_neurons_all_clustering(frame_i,:) = select_neurons_all_clustering{sti_i, session_i, trial_i, frame_i};
            end
            imagesc(curr_select_neurons_all_degree)
            plot([param.fs * param_for_edge.show_sti_before, param.fs * param_for_edge.show_sti_before], [0, size(curr_select_neurons_all_degree,1)], '--r', 'LineWidth', 0.5);
            plot([param.fs * (param_for_edge.show_sti_before + param_for_edge.sti_during_time), param.fs * (param_for_edge.show_sti_before + param_for_edge.sti_during_time)], [0, size(curr_select_neurons_all_degree,1)], '--r', 'LineWidth', 0.5);
            colorbar           
        end

        title(sprintf('Sti %d, Session %d', sti_i, session_i), 'FontSize', 12);
        xlabel('Seconds', 'FontSize', 10);  
        xticklabels = xticks / param.fs - param.show_sti_before; % 将刻度转换为秒数        
        set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);        
        xlabel('Frames', 'FontSize', 10);
        ylabel('neurons', 'FontSize', 10);
        grid off; 
        set(gca, 'FontSize', 10); 

        hold off; 
    end
end
savefig([result_for_response_to_stis_path sprintf('degree_sequence_heatmap')]);
exportgraphics(gcf, [result_for_response_to_stis_path sprintf('degree_sequence_heatmap.png')], 'Resolution', 300);

% catch
% end


