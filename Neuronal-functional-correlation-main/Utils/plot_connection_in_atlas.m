function plot_connection_in_atlas(corr_matrix_cell,select_neurons_C,select_neurons_id,atlas,xlabel_cell,brain_region,neuron_position_x_mat,neuron_position_y_mat,video_name,color)

top_projection_rot = ~rot90(atlas.top_projection,3);
clean_cortex_outline_rot = ~rot90(atlas.clean_cortex_outline,3);
for i =1:length(atlas.parents)
    ids_mat(i)=str2num(atlas.ids{i});
end
for i =2:length(atlas.parents)
    parents_mat(i)=str2num(atlas.parents{i});
end
unique_brain_region_id = [];
region_position_x = zeros(1,length(xlabel_cell));
region_position_y = zeros(1,length(xlabel_cell));
for i = 1 : length(xlabel_cell)
    try
    curr_region = xlabel_cell{i};
    buf = strcmp(brain_region, curr_region);
    C_array{i}= C_mat(find(buf==1), :);
    neuron_num(i) = length(find(buf));
    curr_id = cellfun(@(x) ~isempty(strfind(x,curr_region)),atlas.names);
    unique_brain_region_id(i) = atlas.ids(find(curr_id==1)); 
    neuron_position_x{i} = neuron_position_x_mat(find(buf));
    neuron_position_y{i} = neuron_position_y_mat(find(buf));
    curr_unique_braion_region_id = atlas.ids(find(curr_id==1));
    unique_region_id(i) = find(atlas.unique_top_projection==str2num(curr_unique_braion_region_id{1}));
    unique_region_x(i) = atlas.unique_top_projection_x(unique_region_id(i));
    unique_region_y(i) = atlas.unique_top_projection_y(unique_region_id(i));
    catch
    end
end

mov=VideoWriter(video_name);
mov.FrameRate = 1;
open(mov);
step = 1;
fig = figure;
set(fig,'DoubleBuffer','on');
set(gca,'NextPlot','replace','Visible','off');
for framei = 1:step:size(corr_matrix_cell,2)
    F = corr_matrix_cell{framei};

    %% figure connect in atlas 
    clean_cortex_outline_rot = ~rot90(atlas.clean_cortex_outline,3);

    % dilate atlas
    SE = strel('disk',1); 
    dilated_clean_cortex_outline_rot = imdilate(~clean_cortex_outline_rot,SE);
    
    % adjust color
    dilated_clean_cortex_outline_rot_color = dilated_clean_cortex_outline_rot*0.88;
    dilated_clean_cortex_outline_rot_color(dilated_clean_cortex_outline_rot_color == 0) = 1;

    imshow(dilated_clean_cortex_outline_rot_color);
    hold on;

    select_neurons_CNorm = (select_neurons_C - min(select_neurons_C(:))) / (max(select_neurons_C(:)) - min(select_neurons_C(:)));
    % >>>>>> figure neurons connection
    for kk = 1:length(select_neurons_id)
        hScatter =scatter(neuron_position_x_mat(select_neurons_id(kk)),neuron_position_y_mat(select_neurons_id(kk)),select_neurons_CNorm(kk,framei)*80, color, 'filled');
        hScatter.MarkerFaceAlpha = 0.6; % 调整填充颜色的透明度为50%
        hScatter.MarkerEdgeAlpha = 0.3; % 调整边缘颜色的透明度为50%
        hold on; 
    end

    for i = 1:size(F,1)
        for j = 1:size(F,2)
            aa = neuron_position_x_mat(select_neurons_id(i));
            bb = neuron_position_y_mat(select_neurons_id(i));
            cc = neuron_position_x_mat(select_neurons_id(j));
            dd = neuron_position_y_mat(select_neurons_id(j));
            if F(i,j)>0.9

            x = [aa, cc]; 
            y = [bb, dd]; 
            z = zeros(size(x));     
            patch('XData',x,'YData',y,'ZData',z,'FaceColor','none','EdgeColor',color,'LineWidth',F(i,j),'EdgeAlpha',0.04);
            
            hold on;
            end
        end
    end  
    F_fig = getframe(gca);
    F_fig.cdata = F_fig.cdata(:,:,1); %// 调整视频宽高：H为行数(高)，W为列数(宽)
end
close(mov);   
hold on;
