function [musk] = select_edge_transfer_to_musk(curr_select_edge_id,ori_edge,corr_matrix_point_num)
corr_matrix_id = 1:size(ori_edge,1);
% get triangle musk
triangle_musk = zeros(corr_matrix_point_num,corr_matrix_point_num);
index = 1;
for R_i = 2:corr_matrix_point_num
    for R_j = 1:R_i-1
        triangle_musk(R_i,R_j) = corr_matrix_id(index);
        index = index + 1;
    end
end

musk = zeros(corr_matrix_point_num,corr_matrix_point_num);
for musk_i = 1:size(curr_select_edge_id,1)
    curr_musk_pixel_i = curr_select_edge_id(musk_i);
    [curr_row, curr_column] = find(triangle_musk==curr_musk_pixel_i);
    musk(curr_row, curr_column) = 1;
end
