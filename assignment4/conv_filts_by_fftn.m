%%  convolve filters by fftn
% The temporal filters f1 (fast) and f2 (slow) are convolved with the
% spatial filters oddFilt and evenFilt to produce oddFast, oddSlow,
% evenFast, and evenSlow at each time-point. This algorithm uses fftn to
% compute 3-dimensional fft's for each filter and then computes ifftn of
% the element-wise dot-product of the fft's of spatial and temporal
% filters. The resulting array is 3-dimensional with the third-axis being
% time.
%%
function [oddFast, oddSlow, evenFast, evenSlow] = ...
    conv_filts_by_fftn(f1, f2, oddFilt, evenFilt)

    oddFast = ifftn(fftn(f1) .* fftn(oddFilt));
    oddSlow = ifftn(fftn(f2) .* fftn(oddFilt));
    evenSlow = ifftn(fftn(f2) .* fftn(evenFilt));
    evenFast = ifftn(fftn(f1) .* fftn(evenFilt));
end