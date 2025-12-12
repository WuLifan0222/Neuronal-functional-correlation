function deconve_trace(trace)
neurons_align_deconv = rush_spike_deconv(trace);
figure('position', [1, 1, 1500, 200]); 
imagesc(-zscore(neurons_align_deconv,0, 2));
colormap('gray'); caxis([-2, 0]); axis off;
