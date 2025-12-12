function C_deconv = rush_spike_deconv(C)
C_deconv = zeros(size(C));
for i = 1 : size(C, 1)
    y = C(i, :);
     [c, s, kernel, iter] = deconvCa(y);
     C_deconv(i, :) = c;
end
