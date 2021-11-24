clear; close all; clc;

%%
g = gabor([50 50], [30 330], 'SpatialFrequencyBandwidth', [0.1 0.5]);
%subplot(2, 2, 1);
for p = 1:length(g)
    figure();

    %subplot(2, 2, p)
    imshow(real(g(p).SpatialKernel), []);
end