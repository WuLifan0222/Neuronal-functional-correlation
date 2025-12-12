function [deact_select_bool,deact_response_trials_bool,deact_resp_nums] = estimate_deact_response_by_signrank(data,before,during,after,base_thre,res_nums_thre)
trials_num = size(data,1);
deact_response_trials_bool = zeros(1,trials_num); %
before_stim = data(:,1:before);
after_stim = data(:,before+1:before+during+after);
for triali = 1:trials_num
    if (signrank(after_stim(triali,:),mean(before_stim(triali,:)),'tail','left')<0.05)&&(signrank(after_stim(triali,:),base_thre,'tail','left')<0.05) 
        deact_response_trials_bool(1,triali)=1;
    else
    end
end
deact_resp_nums = sum(deact_response_trials_bool);
deact_select_bool = (deact_resp_nums>=res_nums_thre);
