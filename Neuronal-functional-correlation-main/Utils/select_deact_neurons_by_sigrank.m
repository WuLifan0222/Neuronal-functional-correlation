function [select_deact_neurons_id,select_deact_neurons_C,select_deact_bool_for_all_neurons,deact_response_condition] = select_deact_neurons_by_sigrank(C,stimuli,param)
neuronn = size(C,1);
for neuroni = 1:neuronn      
    % >>>>>>>>> get trials matrix cell for each sti
    [trials_cell,before,during,after,sti_num] = get_trials_matrix_by_sti(C(neuroni,:),stimuli,param);
    % >>>>>>>>> select neurons by signrank
    res_nums_thre = param.res_nums_thre;
    base_thre = mean(C(neuroni,:));
    for i = 1:sti_num
        trials_matrix = trials_cell{1,i};
        [select_bool,response_trials_bool,resp_nums] = estimate_deact_response_by_signrank(trials_matrix,before,during,after,base_thre,res_nums_thre);
        deact_response_condition{i,neuroni}.response_trials_bool = response_trials_bool;
        deact_response_condition{i,neuroni}.resp_nums = resp_nums;
        select_deact_bool_for_all_neurons(i,neuroni) = select_bool;
    end
end
select_deact_neurons_id = find(sum(select_deact_bool_for_all_neurons,1)>=1);
% select_neurons_id = find((sum(select_bool_for_all_neurons,1)>=1) & (sum(select_bool_for_all_neurons,1)<=2)); % adjust here 3
% select_neurons_id = find(sum(select_bool_for_all_neurons,1)==2); % adjust here 3
select_deact_neurons_C = C(select_deact_neurons_id,:);
