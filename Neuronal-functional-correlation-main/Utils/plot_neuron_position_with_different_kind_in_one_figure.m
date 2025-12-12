function plot_neuron_position_with_different_kind_in_one_figure(atlas_bg, resize_ratio, neuron_position_x, neuron_position_y, color, sz)
% background image
atals_mask = atlas_bg;
atals_mask(atals_mask == 0) = 1;
atals_mask(isnan(atals_mask)) = 0;
imagesc(atlas_bg, 'Alphadata', atals_mask); colormap(gray(3)) % note for imagesc, there is no need for change the order

axis equal, axis off
for j = 1 : length(neuron_position_x)
    % processing due to manner of pcolor
    hold on, scatter((neuron_position_x(j) * resize_ratio), (neuron_position_y(j) * resize_ratio), ...
                    sz, color, 'filled', 'MarkerFaceAlpha',.8,'MarkerEdgeAlpha',.8)               
end
