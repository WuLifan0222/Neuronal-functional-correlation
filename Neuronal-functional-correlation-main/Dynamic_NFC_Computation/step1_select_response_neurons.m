clc; clear; close all; set(0,'defaultfigurecolor','w');

%% analysis sti response --WLF 20230811
addpath(genpath('D:\PPT_gather\Neuronal functional connectivity\Utils'));

%% >>>>>>>>> load color
color_scheme_npg = get_color();

%% >>>>>>>>> Load experiment 
folder_path = ['D:\PPT_gather\Neuronal functional connectivity\dataset\'];
load([folder_path 'demo_neurons_data.mat']);
C = C_mat;

%% stimuli
% >>>>>>>>> find different stis
load([folder_path 'demo_visual_stimuli_with_label.mat']);
stimuli = load([folder_path 'demo_visual_stimuli_with_label.mat']);

% >>>>>>>>>>>> preprocess for different lable
buff = min(stimuli.stimuili_label_ind)-1;
stimuli.stimuili_label_ind = stimuli.stimuili_label_ind-buff;
stimuli.stimuili_label = stimuli.stimuili_label_ind';
stimuli.stimuli_array_with_label = stimuli.stimuli_array_with_label-buff;
stimuli.stimuli_array_with_label(stimuli.stimuli_array_with_label < 0) = 0;
tabulated_for_sti = tabulate(stimuli.stimuili_label_ind);

% >>>>>> process for null sti num
sti_num = size(tabulated_for_sti,1);
sti_num = find(tabulated_for_sti(:,2)~=0);
% trial_nums = min(tabulated_for_sti(:,2));
trial_nums = 8;  %%%  1-8th for select 9-20 for classification task
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

%% >>>>>>>>> atlas
atlas_bg = importdata('D:\PPT_gather\Neuronal functional connectivity\Utils\atlas\cortical_out_line_resize_5_30.mat');
atlas = importdata('D:\PPT_gather\Neuronal functional connectivity\Utils\atlas\atlas_top_projection_30.mat');
resize_ratio = 5;
% atlas erode for show 
matrix_nan_to_zero = atlas_bg;
matrix_nan_to_zero(matrix_nan_to_zero == 0) = 1;
matrix_nan_to_zero(isnan(matrix_nan_to_zero)) = 0;
se = strel('disk', 8); 
eroded_matrix = imerode(matrix_nan_to_zero, se);
eroded_matrix_with_nan = eroded_matrix;
eroded_matrix_with_nan(eroded_matrix_with_nan == 0) = NaN;
eroded_matrix_with_nan(eroded_matrix_with_nan == 1) = 0;
atlas_bg = eroded_matrix_with_nan;


%% >>>>>>>>> set experiment parameters
param.fs = 4;
param.sti_during_time = 3; %change here according different sti
param.show_sti_before = 3; %change here according different sti
param.show_sti_after = 7; %change here according different sti
param.delay_after = 0;
param.res_nums_thre = 8; % change here select neurons by response nums -- response reliability:(8/8=100%)
param.show_fig_for_all_neurons = 0; % 1:show each neuron trace/trials imagesc
param.show_fig_for_select_neurons = 0; % 1:show each neuron trace/trials imagesc
param.show_fig_for_select_good_neurons = 0; % 1:show selected osi>0.8 neurons
delay_win = param.sti_during_time + param.show_sti_before +param.show_sti_after;
% auto get more param
[trials_cell, before, during, after, sti_labels, sti_num, trials_num_matrix] = get_multi_neurons_trials_matrix_by_sti(C(1,:), stimuli, param);
param.sti_labels = sti_labels;
param.sti_num = sti_num;
param.trials_num = min(trials_num_matrix);

% easy select by trials intensity & sort by similarity
param.trial_select_thresh = 1.7; % during intensity > before intensity * param.trial_select_thresh
param.response_ratio = 0.8;
[C_select_by_trials_ind,C_select_by_trials,C_trials,neuron_trial_similar,select_response_by_trial_for_each_sti_condition] = select_by_trials(C,param,stimuli);

%% >>>>>>>>> select neuron for all brain region
figure_path = [folder_path sprintf('result_for_thre%d\',param.res_nums_thre)];
result_for_response_to_stis_path = [figure_path '\'];
mkdir(result_for_response_to_stis_path);

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  act 
C_select_for_multi_sti_cell_ind = {};
C_select_for_multi_sti_cell = {};
[select_neurons_id,select_neurons_C,select_bool_for_all_neurons,response_condition] = select_neurons_by_sigrank(C,stimuli,param);


% >>>>>> show each selected neurons
if param.show_fig_for_select_neurons == 1
for select_neuron_i = 1:size(select_neurons_id,2)
    figure(Position=[100,100,700,800])
    sgtitle(sprintf('response of neuron %d',select_neurons_id(select_neuron_i)));
    for sti_i = 1:sti_num
        single_trial =C_trials{select_neurons_id(select_neuron_i)}{sti_i}';
        neuron_single_trial_average = mean(single_trial);
        subplot(2,sti_num,sti_i);
        imagesc(single_trial);hold on;
        plot([param.fs*(param.show_sti_before),param.fs*(param.show_sti_before)],[0,size(single_trial,1)],'--r','linewidth',0.5);hold on;
        plot([param.fs*(param.show_sti_before+param.sti_during_time),param.fs*(param.show_sti_before+param.sti_during_time)],[0,size(single_trial,1)],'--r','linewidth',0.5);hold on;            
        ylabel('trials #');
        subplot(2,sti_num,sti_i+sti_num);
        hold on;plot(single_trial','Color', [0.7 0.7 0.7]);
        plot(neuron_single_trial_average,'Color',color_scheme_npg(1,:),'LineWidth',2);ylabel('intensity'); hold on;
        plot([param.fs*(param.show_sti_before),param.fs*(param.show_sti_before)],[-max(max(single_trial)),2*max(max(single_trial))],'--r','linewidth',0.5);hold on;
        plot([param.fs*(param.show_sti_before+param.sti_during_time),param.fs*(param.show_sti_before+param.sti_during_time)],[-max(max(single_trial)),2*max(max(single_trial))],'--r','linewidth',0.5);hold on;      
        xlabel('frame (4hz)');          
    end
    savefig([result_for_response_to_stis_path sprintf('select_neuron%d',select_neurons_id(select_neuron_i)) ]);
    exportgraphics(gcf,[result_for_response_to_stis_path sprintf('select_neuron%d',select_neurons_id(select_neuron_i))  '.png'],'Resolution',300)
    select_neuron_i
    close all
end
end

select_response_by_sigrank_for_each_sti_condition = select_bool_for_all_neurons;
response_to_stis_neurons = sum(select_bool_for_all_neurons,2);
all_neuron_response_times = sum(select_bool_for_all_neurons,1);


% >>>>>>>>>>>> neurons response to different stis
for i = 1:sti_num
    response_neuron_id = find(select_bool_for_all_neurons(i,:)==1);
    % >>>>>> plot trace
    response_neuron_C = C(response_neuron_id,:);

    if param.show_fig_for_select_neurons == 1
    % >>>>>> show each selected neurons
    for select_neuron_i = 1:size(response_neuron_C,1)
        figure(Position=[100,100,700,800])
        sgtitle(sprintf('response of neuron %d',response_neuron_id(select_neuron_i)));
        for sti_i = 1:sti_num
            single_trial =C_trials{response_neuron_id(select_neuron_i)}{sti_i}';
            neuron_single_trial_average = mean(single_trial);
            subplot(2,sti_num,sti_i);
            imagesc(single_trial);hold on;
            plot([param.fs*(param.show_sti_before),param.fs*(param.show_sti_before)],[0,size(single_trial,1)],'--r','linewidth',0.5);hold on;
            plot([param.fs*(param.show_sti_before+param.sti_during_time),param.fs*(param.show_sti_before+param.sti_during_time)],[0,size(single_trial,1)],'--r','linewidth',0.5);hold on;            
            ylabel('trials #');
            subplot(2,sti_num,sti_i+sti_num);
            hold on;plot(single_trial','Color', [0.7 0.7 0.7]);
            plot(neuron_single_trial_average,'Color',color_scheme_npg(1,:),'LineWidth',2);ylabel('intensity'); hold on;
            plot([param.fs*(param.show_sti_before),param.fs*(param.show_sti_before)],[-max(max(single_trial)),2*max(max(single_trial))],'--r','linewidth',0.5);hold on;
            plot([param.fs*(param.show_sti_before+param.sti_during_time),param.fs*(param.show_sti_before+param.sti_during_time)],[-max(max(single_trial)),2*max(max(single_trial))],'--r','linewidth',0.5);hold on;      
            xlabel('frame (4hz)');          
        end
        savefig([result_for_response_to_stis_path sprintf('select_neuron%d',response_neuron_id(select_neuron_i)) ]);
        exportgraphics(gcf,[result_for_response_to_stis_path sprintf('select_neuron%d',response_neuron_id(select_neuron_i))  '.png'],'Resolution',300)
        select_neuron_i
        close all
    end
    end
end


% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  deact
C_select_deact_for_multi_sti_cell_ind = {};
C_select_deact_for_multi_sti_cell = {};
[select_neurons_id,select_neurons_C,select_bool_for_all_neurons,response_condition] = select_neurons_by_sigrank(C,stimuli,param);
[select_deact_neurons_id,select_deact_neurons_C,select_deact_bool_for_all_neurons,deact_response_condition] = select_deact_neurons_by_sigrank(C,stimuli,param);

if param.show_fig_for_select_neurons == 1
% >>>>>> show each selected neurons
for select_neuron_i = 1:size(select_deact_neurons_id,2)
    figure(Position=[100,100,700,800])
    sgtitle(sprintf('response of neuron %d',select_deact_neurons_id(select_neuron_i)));
    for sti_i = 1:sti_num
        single_trial =C_trials{select_deact_neurons_id(select_neuron_i)}{sti_i}';
        neuron_single_trial_average = mean(single_trial);
        subplot(2,sti_num,sti_i);
        imagesc(single_trial);hold on;
        plot([param.fs*(param.show_sti_before),param.fs*(param.show_sti_before)],[0,size(single_trial,1)],'--r','linewidth',0.5);hold on;
        plot([param.fs*(param.show_sti_before+param.sti_during_time),param.fs*(param.show_sti_before+param.sti_during_time)],[0,size(single_trial,1)],'--r','linewidth',0.5);hold on;            
        ylabel('trials #');
        subplot(2,sti_num,sti_i+sti_num);
        hold on;plot(single_trial','Color', [0.7 0.7 0.7]);
        plot(neuron_single_trial_average,'Color',color_scheme_npg(1,:),'LineWidth',2);ylabel('intensity'); hold on;
        plot([param.fs*(param.show_sti_before),param.fs*(param.show_sti_before)],[-max(max(single_trial)),2*max(max(single_trial))],'--r','linewidth',0.5);hold on;
        plot([param.fs*(param.show_sti_before+param.sti_during_time),param.fs*(param.show_sti_before+param.sti_during_time)],[-max(max(single_trial)),2*max(max(single_trial))],'--r','linewidth',0.5);hold on;      
        xlabel('frame (4hz)');          
    end
    savefig([result_for_response_to_stis_path sprintf('select_neuron%d',select_deact_neurons_id(select_neuron_i)) ]);
    exportgraphics(gcf,[result_for_response_to_stis_path sprintf('select_neuron%d',select_deact_neurons_id(select_neuron_i))  '.png'],'Resolution',300)
    select_neuron_i
    close all
end
end

select_deact_response_by_sigrank_for_each_sti_condition = select_deact_bool_for_all_neurons;
deact_response_to_stis_neurons = sum(select_deact_bool_for_all_neurons,2);
all_neuron_deact_response_times = sum(select_deact_bool_for_all_neurons,1);

% >>>>>>>>>>>> neurons response to different stis
for i = 1:sti_num
    deact_response_neuron_id = find(select_deact_bool_for_all_neurons(i,:)==1);
    % >>>>>> plot trace
    deact_response_neuron_C = C(deact_response_neuron_id,:);
    C_select_deact_for_multi_sti_cell_ind{i} = deact_response_neuron_id;
    C_select_deact_for_multi_sti_cell{i} = deact_response_neuron_C;

    % >>>>>> show each deact selected neurons
    if param.show_fig_for_select_neurons == 1
    for deact_select_neuron_i = 1:size(deact_response_neuron_C,1)
        figure(Position=[100,100,700,800])
        sgtitle(sprintf('deact response of neuron %d',deact_response_neuron_id(deact_select_neuron_i)));
        for sti_i = 1:sti_num
            deact_single_trial =C_trials{deact_response_neuron_id(deact_select_neuron_i)}{sti_i}';
            deact_neuron_single_trial_average = mean(deact_single_trial);
            subplot(2,sti_num,sti_i);
            imagesc(deact_single_trial);hold on;
            plot([param.fs*(param.show_sti_before),param.fs*(param.show_sti_before)],[0,size(deact_single_trial,1)],'--r','linewidth',0.5);hold on;
            plot([param.fs*(param.show_sti_before+param.sti_during_time),param.fs*(param.show_sti_before+param.sti_during_time)],[0,size(deact_single_trial,1)],'--r','linewidth',0.5);hold on;            
            ylabel('trials #');
            subplot(2,sti_num,sti_i+sti_num);
            hold on;plot(deact_single_trial','Color', [0.7 0.7 0.7]);
            plot(deact_neuron_single_trial_average,'Color',color_scheme_npg(1,:),'LineWidth',2);ylabel('intensity'); hold on;
            plot([param.fs*(param.show_sti_before),param.fs*(param.show_sti_before)],[-max(max(deact_single_trial)),2*max(max(deact_single_trial))],'--r','linewidth',0.5);hold on;
            plot([param.fs*(param.show_sti_before+param.sti_during_time),param.fs*(param.show_sti_before+param.sti_during_time)],[-max(max(deact_single_trial)),2*max(max(deact_single_trial))],'--r','linewidth',0.5);hold on;      
            xlabel('frame (4hz)');          
        end
        savefig([result_for_response_to_stis_path sprintf('deact_select_neuron%d',deact_response_neuron_id(deact_select_neuron_i)) ]);
        exportgraphics(gcf,[result_for_response_to_stis_path sprintf('deact_select_neuron%d',deact_response_neuron_id(deact_select_neuron_i)) '.png'],'Resolution',300)
    end
    end
end


%% show act deact together
% >>>>> figure in one fig
figure('position', [100, 100, 600, 500]), 
axes('position', [0, 0, 1, 1])
for i = 1:sti_num
    act_response_neuron_id = find(select_bool_for_all_neurons(i,:)==1);
    deact_response_neuron_id = find(select_deact_bool_for_all_neurons(i,:)==1);
    % plot neuron position distribution
    plot_neuron_position_with_different_kind_in_one_figure(atlas_bg, resize_ratio, neuron_position_x_mat(act_response_neuron_id), ...
        neuron_position_y_mat(act_response_neuron_id), color_scheme_npg(1, :), 30)
    hold on
    plot_neuron_position_with_different_kind_in_one_figure(atlas_bg, resize_ratio, neuron_position_x_mat(deact_response_neuron_id), ...
        neuron_position_y_mat(deact_response_neuron_id), color_scheme_npg(2, :), 30)
    hold on
end

figure
% plot neuron position distribution
plot_neuron_position_with_different_kind_in_one_figure(atlas_bg, resize_ratio, neuron_position_x_mat(select_neurons_id), ...
    neuron_position_y_mat(select_neurons_id), color_scheme_npg(1, :), 30)
hold on
plot_neuron_position_with_different_kind_in_one_figure(atlas_bg, resize_ratio, neuron_position_x_mat(select_deact_neurons_id), ...
    neuron_position_y_mat(select_deact_neurons_id), color_scheme_npg(2, :), 30)
hold on

save([result_for_response_to_stis_path 'select_neurons.mat'], ...
    'C','neuron_position_x_mat','neuron_position_y_mat','C_trials','param','stimuli',...
    'select_neurons_id','select_neurons_C','C_select_for_multi_sti_cell_ind', 'C_select_for_multi_sti_cell',...
    'select_bool_for_all_neurons','response_condition','select_response_by_sigrank_for_each_sti_condition',...
    'select_deact_neurons_id','select_deact_neurons_C','C_select_deact_for_multi_sti_cell_ind','C_select_deact_for_multi_sti_cell', ...
    'select_deact_bool_for_all_neurons','deact_response_condition','select_deact_response_by_sigrank_for_each_sti_condition','-v7.3');

save([folder_path 'select_neurons.mat'], ...
    'C','neuron_position_x_mat','neuron_position_y_mat','C_trials','param','stimuli',...
    'select_neurons_id','select_neurons_C','C_select_for_multi_sti_cell_ind', 'C_select_for_multi_sti_cell',...
    'select_bool_for_all_neurons','response_condition','select_response_by_sigrank_for_each_sti_condition',...
    'select_deact_neurons_id','select_deact_neurons_C','C_select_deact_for_multi_sti_cell_ind','C_select_deact_for_multi_sti_cell', ...
    'select_deact_bool_for_all_neurons','deact_response_condition','select_deact_response_by_sigrank_for_each_sti_condition','-v7.3');
