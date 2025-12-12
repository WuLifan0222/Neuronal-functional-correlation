function plot_multi_trace_with_one_sti_labels(trace, stimuli, color_scheme_npg,stii)
start_edge = stimuli.start_edge;
end_edge = stimuli.end_edge;
stimuili_label_ind = stimuli.stimuili_label_ind;
stimuili_label = stimuli.stimuili_label;
temporal_trace_render(trace)  

% add stimuli as background shadow
for i = 1 : length(start_edge) 
    if (stimuili_label(i)==stii)
        hold on, rectangle('Position', [start_edge(i),-1, end_edge(i) - start_edge(i),size(trace, 1) * 3 + 2],...
            'facecolor', [color_scheme_npg(stimuili_label_ind(i), :), 0.4], 'edgecolor', [color_scheme_npg(stimuili_label_ind(i), :), 0.4], ...
            'linewidth', 0.2);hold on 
        try
        txt = num2str(stimuili_label_ind(1,i));
        catch
        txt = num2str(stimuili_label(1,i));
        end
        text(start_edge(i)-1,-3,txt);hold on 
    else 
    end
end
