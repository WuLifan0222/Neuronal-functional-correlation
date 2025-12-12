function mean_dist = compute_matrixs_similarity(matrix_cell)
n = size(matrix_cell,2);  % 获取矩阵个数
total_dist = 0; 
count = 0;
% 遍历所有矩阵对
for i = 1:n-1
  for j = i+1:n
    % 计算两个矩阵之间的欧几里得距离
    matrix_i = matrix_cell{i}; matrix_i(isnan(matrix_i)) = 0;
    matrix_j = matrix_cell{j}; matrix_j(isnan(matrix_j)) = 0;
    dist = norm(matrix_i - matrix_j, 'fro'); 
    % 累加到总距离
    total_dist = total_dist + dist; 
    count = count + 1;
  end
end
mean_dist = total_dist/count;
