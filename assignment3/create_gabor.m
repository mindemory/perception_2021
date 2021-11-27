clear; close all; clc;

%%
g = gabor([10 10], [30 330], 'SpatialFrequencyBandwidth', [0.3 0.5]);
%subplot(2, 2, 1);
figure();
for p = 1:length(g)
    %figure();

    subplot(2, 2, p)
    I = real(g(p).SpatialKernel);
    imshow(I, [min(I(:)) * 2, max(I(:))]);
end

