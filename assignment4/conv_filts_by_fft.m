%%  convolve filters by fft
% The temporal filters f1 (fast) and f2 (slow) are convolved with the
% spatial filters oddFilt and evenFilt to produce oddFast, oddSlow,
% evenFast, and evenSlow at each time-point. This algorithm uses fft2 to
% compute 2-dimensional fft's for each filter and then computes ifft2 of
% the element-wise dot-product of the fft's of spatial and temporal
% filters at each time-point.
%%
function [oddFast, oddSlow, evenFast, evenSlow] = ...
    conv_filts_by_fft(f1, f2, oddFilt, evenFilt)
    [lx, ly, lt] = size(f1);
    oddFast = zeros(lx, ly, lt);
    oddSlow = zeros(lx, ly, lt);
    evenSlow = zeros(lx, ly, lt);
    evenFast = zeros(lx, ly, lt);
    
    for tt = 1:lt
        oddFast(:, :, tt) = ifft2(fft2(f1(:, :, tt)) .* fft2(oddFilt));
        oddSlow(:, :, tt) = ifft2(fft2(f2(:, :, tt)) .* fft2(oddFilt));
        evenSlow(:, :, tt) = ifft2(fft2(f2(:, :, tt)) .* fft2(evenFilt));
        evenFast(:, :, tt) = ifft2(fft2(f1(:, :, tt)) .* fft2(evenFilt));
    end
end
