function [similarity_index,mean_similarity,degree_sequence,closeness_sequence,betweenness_sequence,average_path_length_sequence,clustering_sequence,...
all_degree,all_closeness,all_betweenness,all_average_path_length,all_clustering] = ...
graph_network_analysis(C_mat,select_neurons_C_sort_by_region_edge_cell,stimuli,stimuli_for_edge,param)

[corr_trials_cell,before,during,after,sti_num] = get_corr_trials_matrix_by_sti(select_neurons_C_sort_by_region_edge_cell,stimuli,param);
C_matNorm = (C_mat - min(C_mat(:))) / (max(C_mat(:)) - min(C_mat(:)));
[C_matNorm_cell,before,during,after,sti_num] = get_trials_matrix_multi_trace_by_sti((C_matNorm+0.0000001),stimuli,param);

each_session_trials_num = param.trials_num/param.data_session_n;

% compute some index 
all_degree = {};
all_closeness = {};
all_betweenness = {};
all_average_path_length = {};
all_clustering = {};

degree_sequence = [];
closeness_sequence = [];
betweenness_sequence = [];
average_path_length_sequence = [];
clustering_sequence = [];

for sti_i = 1:param.sti_num
    for data_sessioni = 1:param.data_session_n
        curr_matrix_all_trial = corr_trials_cell{sti_i};
        n = 0;
        for trial_i = each_session_trials_num*(data_sessioni-1)+1:each_session_trials_num*data_sessioni
            n = n+1;
            curr_matrix_to_compute_all{n} = curr_matrix_all_trial{trial_i,24}; %% change here 24 means the24th frame    
            for frame_i = 1: size(curr_matrix_all_trial,2)
                curr_matrix_to_compute_allframe_all{n,frame_i} = curr_matrix_all_trial{trial_i,frame_i};
            end
        end
        % index1 : trial musk corr matrix similarity (for 24th frames)
        similarity_index(data_sessioni,sti_i) = compute_matrixs_similarity(curr_matrix_to_compute_all);

        for one_session_trial_i = 1:size(curr_matrix_to_compute_allframe_all,1)
            for frame_i = 1: size(curr_matrix_to_compute_allframe_all,2)
                curr_matrix = curr_matrix_to_compute_allframe_all{one_session_trial_i,frame_i};
                curr_matrix(isnan(curr_matrix) | isinf(curr_matrix)) = 0;
                curr_matrix = curr_matrix > param.edeg_threshold;
                brain_graph = graph(curr_matrix);

                % 计算连通组件数量
                components = conncomp(brain_graph); % 找到连接组件
                num_components = max(components); % 计算组件数量

                % 计算平均路径长度（如果是连通的）
                if num_components == 1
                    meanPathLength = mean(distances(brain_graph));
                else
                    % disp('网络不是连通的，不能计算平均路径长度');
                    meanPathLength = NaN; 
                end

                % 计算聚集系数
                clusteringCoefficients = clustering_coefficients(brain_graph);
                meanClusteringCoefficient = mean(clusteringCoefficients);

                % 计算图论特征
                degree = centrality(brain_graph, 'degree');
                closeness = centrality(brain_graph, 'closeness');
                betweenness = centrality(brain_graph, 'betweenness');

                all_average_path_length{sti_i,data_sessioni,one_session_trial_i,frame_i} =meanPathLength;
                all_clustering{sti_i,data_sessioni,one_session_trial_i,frame_i} = clusteringCoefficients;
                all_degree{sti_i,data_sessioni,one_session_trial_i,frame_i} = degree;
                all_closeness{sti_i,data_sessioni,one_session_trial_i,frame_i} = closeness;
                all_betweenness{sti_i,data_sessioni,one_session_trial_i,frame_i} = betweenness;             

                average_path_length_sequence(sti_i,data_sessioni,one_session_trial_i,frame_i) = mean(meanPathLength);
                clustering_sequence(sti_i,data_sessioni,one_session_trial_i,frame_i) = meanClusteringCoefficient;
                degree_sequence(sti_i,data_sessioni,one_session_trial_i,frame_i) = mean(degree);
                closeness_sequence(sti_i,data_sessioni,one_session_trial_i,frame_i) = mean(closeness);
                betweenness_sequence(sti_i,data_sessioni,one_session_trial_i,frame_i) = mean(betweenness);
            end
        end
    end   
end
mean_similarity = mean(similarity_index);
    
