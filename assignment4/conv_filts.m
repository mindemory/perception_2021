%%  convolve filters
% The temporal filters f1 (fast) and f2 (slow) are convolved with the
% spatial filters oddFilt and evenFilt to produce oddFast, oddSlow,
% evenFast, and evenSlow at each time-point. This algorithm uses conv2 to
% compute 2-dimensional convolutions at each time-point. 
%%
function [oddFast, oddSlow, evenFast, evenSlow] = ...
    conv_filts(f1, f2, oddFilt, evenFilt)
    [lx, ly, lt] = size(f1);
    oddFast = zeros(lx, ly, lt);
    oddSlow = zeros(lx, ly, lt);
    evenSlow = zeros(lx, ly, lt);
    evenFast = zeros(lx, ly, lt);
    
    for tt = 1:lt
        oddFast(:, :, tt) = conv2(f1(:, :, tt), oddFilt, 'same');
        oddSlow(:, :, tt) = conv2(f2(:, :, tt), oddFilt, 'same');
        evenSlow(:, :, tt) = conv2(f2(:, :, tt), evenFilt, 'same');
        evenFast(:, :, tt) = conv2(f1(:, :, tt), evenFilt, 'same');
    end
end