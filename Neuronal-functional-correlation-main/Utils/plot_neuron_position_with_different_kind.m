function plot_neuron_position_with_different_kind(atlas_bg, resize_ratio, neuron_position_x, neuron_position_y, my_colormap, sz)
figure('position', [100, 100, 600, 500]), 
axes('position', [0, 0, 1, 1])

e1 = pcolor(atlas_bg); ...
shading interp, colormap(gray(1));
axis equal, axis off

for j = 1 : length(neuron_position_x)
    % processing due to manner of pcolor
    hold on, scatter((neuron_position_x(j) * resize_ratio), (neuron_position_y(j) * resize_ratio), ...
                    sz, my_colormap(9, :), 'filled', 'MarkerFaceAlpha',.4,'MarkerEdgeAlpha',.4)               
end
