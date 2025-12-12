function plot_boundedline_for_trials_matrix(trials_matrix,color_scheme_npg_i)
trial_average_of_each_sti = mean(trials_matrix,1);
mm = [];
err = [];
C_stii_mean = mean(trials_matrix,2);
mm=double(mean(trials_matrix,1));
err = double(std(trials_matrix,[],1)/sqrt(size(trials_matrix,1)));            
ax = boundedline(1:1:length(mm),mm,err,'alpha','cmap',color_scheme_npg_i);hold on;
ax.LineWidth = 2; hold on;
