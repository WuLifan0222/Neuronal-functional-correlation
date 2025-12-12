function [avg_matrix_sort,sortedTimes, sortedValues,sortIdx] = sort_trace_by_peakTime(avg_matrix)
% sort peaks
% 初始化两个向量，用于存储峰值及其出现时间
peakValues = zeros(1, size(avg_matrix,2));         % 保存每列的峰值
peakTimes = zeros(1, size(avg_matrix,2));          % 保存每列峰值出现的行索引（时间点）

% 遍历每列，找到峰值及其对应索引（时间点）
for col = 1:size(avg_matrix,2)
    [peakValues(col), peakTimes(col)] = max(avg_matrix(:, col));  % 找到峰值及其索引
end

% 根据峰值出现的时间点进行排序
[sortedTimes, sortIdx] = sort(peakTimes);  % 排序时间从小到大
sortedValues = peakValues(sortIdx);       % 对应的排序后峰值
avg_matrix_sort = avg_matrix(:,sortIdx);  

