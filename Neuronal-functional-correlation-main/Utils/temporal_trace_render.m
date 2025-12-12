function temporal_trace_render(C)
    ts = zscore(C , 0, 2);
    y_shift = 3;
    clip = false;
    sel = 1:size(ts,1);
    nixs = 1:size(ts,1);
    sel_nixs = nixs(sel);

    length_t = sqrt(size(C, 2)) * 10;
    length_n = sqrt(size(C, 1)) * 60;
    t = (0:size(ts,2)-1);

    figure('Position', [100, 100, length_t, length_n ]);
    hold on
    for n_ix = 1:numel(sel_nixs)
        ax = gca();
        ax.ColorOrderIndex = 1;
        loop_ts = ts(sel_nixs(n_ix),:);
        if clip
            loop_ts(loop_ts > 3*y_shift) = y_shift;
            loop_ts(loop_ts < -3*y_shift) = -y_shift;
        end

        if max(loop_ts) - mean(loop_ts) > 2.5 * y_shift       
            plot1 = plot(t, squeeze(loop_ts) + y_shift*(n_ix-1), 'color', 'k');
            plot1.Color(4) = 0.8;
        else
            % with alpha
            plot1 = plot(t, squeeze(loop_ts) + y_shift*(n_ix-1), 'k');
            plot1.Color(4) = 0.8;
        end   
    end
    xlabel('Frame');
    xlim([min(t) max(t)]);
    axis off
    hold off;
    axis tight;
    set(gca,'LooseInset',get(gca,'TightInset'))
end