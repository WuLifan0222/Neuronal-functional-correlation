function deconve_trace_with_label_gray(trace, stimuli, color_scheme_npg)
start_edge = stimuli.start_edge;
end_edge = stimuli.end_edge;
stimuili_label_ind = stimuli.stimuili_label_ind;
stimuili_label = stimuli.stimuili_label;
deconve_trace_gray(trace);  
% add stimuli as background shadow
for i = 1 : length(start_edge)  
    hold on, rectangle('Position', [start_edge(i),-1, end_edge(i) - start_edge(i),size(trace, 1) * 3 + 2],...
        'facecolor', [color_scheme_npg(stimuili_label_ind(i), :), 0.4], 'edgecolor', [color_scheme_npg(stimuili_label_ind(i), :), 0.4], ...
        'linewidth', 0.2);hold on 
    % txt = num2str(stimuili_label_ind(1,i));
    % text(start_edge(i)-1,-3,txt);hold on 
end
