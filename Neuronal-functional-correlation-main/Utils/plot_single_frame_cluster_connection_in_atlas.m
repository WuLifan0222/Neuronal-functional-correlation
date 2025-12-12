function plot_single_frame_cluster_connection_in_atlas(C_mat,F,F_musk,edge_thre,C_all_value,C_value,select_neurons_id,atlas,xlabel_cell,brain_region,neuron_position_x_mat,neuron_position_y_mat,color,text_title)

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
    curr_id = cellfun(@(x) contains(x,curr_region),atlas.names);
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

% figure connect in atlas
clean_cortex_outline_rot = ~(atlas.clean_cortex_outline);
% dilate atlas
SE = strel('disk',1); % 创建一个半径为1的圆形结构元素
dilated_clean_cortex_outline_rot = imdilate(~clean_cortex_outline_rot,SE);

% adjust color
dilated_clean_cortex_outline_rot_color = dilated_clean_cortex_outline_rot*0.88;
dilated_clean_cortex_outline_rot_color(dilated_clean_cortex_outline_rot_color == 0) = 1;

% imshow(dilated_clean_cortex_outline_rot_color(184:end-150,228:end-70));
imshow(dilated_clean_cortex_outline_rot_color);
hold on;

% >>>>>> figure selected neurons connection
for kk = 1:length(select_neurons_id)
% for kk = 1:100:length(select_neurons_id)
    hScatter =scatter(neuron_position_x_mat(select_neurons_id(kk)),neuron_position_y_mat(select_neurons_id(kk)),C_value(kk)*40, color, 'filled');
    hScatter.MarkerFaceAlpha = 0.3; % 调整填充颜色的透明度为50%
    hScatter.MarkerEdgeAlpha = 0.3; % 调整边缘颜色的透明度为50%
    hold on; 
end

for i = 1:size(F,1)
    for j = 1:size(F,2)
        aa = neuron_position_x_mat(select_neurons_id(i));
        bb = neuron_position_y_mat(select_neurons_id(i));
        cc = neuron_position_x_mat(select_neurons_id(j));
        dd = neuron_position_y_mat(select_neurons_id(j));
        % if F(i,j)>0.9
        if F(i,j)>= edge_thre
            x = [aa, cc]; % x坐标
            y = [bb, dd]; % y坐标
            z = zeros(size(x)); % 设置为与x，y同样长度的零矩阵，因为我们绘制的是2D图形      
            patch('XData',x,'YData',y,'ZData',z,'FaceColor','none','EdgeColor',color,'LineWidth',F(i,j),'EdgeAlpha',0.04);              
            hold on;
        end
    end
end  

% >>>>>> figure cluster neurons connection
[row,col] = find(F_musk>edge_thre);  
cluster_id = unique([row,col]);
for kk = 1:length(cluster_id)
    hScatter =scatter(neuron_position_x_mat(select_neurons_id(cluster_id(kk))),neuron_position_y_mat(select_neurons_id(cluster_id(kk))),C_value(cluster_id(kk))*60, color, 'filled');
    hScatter.MarkerFaceAlpha = 0.8; % 调整填充颜色的透明度为50%
    hScatter.MarkerEdgeAlpha = 0.8; % 调整边缘颜色的透明度为50%
    hold on; 
end

for i = 1:size(F_musk,1)
    for j = 1:size(F_musk,2)
        aa = neuron_position_x_mat(select_neurons_id(i));
        bb = neuron_position_y_mat(select_neurons_id(i));
        cc = neuron_position_x_mat(select_neurons_id(j));
        dd = neuron_position_y_mat(select_neurons_id(j));
        if F_musk(i,j)>edge_thre
            x = [aa, cc]; % x坐标
            y = [bb, dd]; % y坐标
            z = zeros(size(x)); % 设置为与x，y同样长度的零矩阵，因为我们绘制的是2D图形      
            patch('XData',x,'YData',y,'ZData',z,'FaceColor','none','EdgeColor',color,'LineWidth',F_musk(i,j)*2,'EdgeAlpha',0.3);  
            hold on;
        end
    end
end 

title(text_title, 'FontWeight', 'bold', 'FontSize', 12, 'FontName', 'Times New Roman');
F_fig = getframe(gca);
F_fig.cdata = F_fig.cdata(:,:,1); %// 调整视频宽高：H为行数(高)，W为列数(宽)
hold on;
