function deconve_trace_with_label(trace, stimuli, color_scheme_npg)
start_edge = stimuli.start_edge;
end_edge = stimuli.end_edge;
stimuili_label_ind = stimuli.stimuili_label_ind;
stimuili_label = stimuli.stimuili_label;
deconve_trace(trace)  
% add stimuli as background shadow
for i = 1 : length(start_edge)  
    hold on, rectangle('Position', [start_edge(i),-size(trace, 1)/10-5, end_edge(i) - start_edge(i), size(trace, 1)/10+7],...
        'facecolor', [color_scheme_npg(stimuili_label_ind(i), :), 1], 'edgecolor', [color_scheme_npg(stimuili_label_ind(i), :), 1], ...
        'linewidth', 1);hold on 
end
