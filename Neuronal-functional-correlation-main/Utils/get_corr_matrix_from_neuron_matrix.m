function [corr_matrix,edge_cell] = get_corr_matrix_from_neuron_matrix(C,corr_param)
corr_win_num_i = 1;
corr_matrix = [];
R_cell = {};
edge_vector_cell = {};

for corr_win_before_num_i = 1:corr_param.corr_win_num   
    edge_vector = [];
    R = corrcoef(C(:,corr_win_num_i:corr_win_num_i+corr_param.corr_win-1)');    
    R_lowerTriangularVector = R(tril(true(size(R)), -1));

    try
    corr_matrix = [corr_matrix,R_lowerTriangularVector];
    catch
    corr_matrix = [corr_matrix,zeros(size(corr_matrix,1),1)+0.000000001];
    end

    edge_cell{corr_win_before_num_i} = R ; 
end
corr_matrix(isnan(corr_matrix)) = 0;
