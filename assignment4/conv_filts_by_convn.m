%%  convolve filters by convn
% The temporal filters f1 (fast) and f2 (slow) are convolved with the
% spatial filters oddFilt and evenFilt to produce oddFast, oddSlow,
% evenFast, and evenSlow at each time-point. This algorithm uses convn to
% compute 3-dimensional convolutions. The resulting array is 3-dimensional 
% with the third-axis being time.
%%
function [oddFast, oddSlow, evenFast, evenSlow] = ...
    conv_filts_by_convn(f1, f2, oddFilt, evenFilt)
    oddFast = convn(f1, oddFilt, 'same');
    oddSlow = convn(f2, oddFilt, 'same');
    evenSlow = convn(f2, evenFilt, 'same');
    evenFast = convn(f1, evenFilt, 'same');
end