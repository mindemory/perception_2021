%%  gabor filter
% Given a spatial resolution x, the standard deviation of the gaussian
% sigma, and the spatial frequency sf, gabor_filter computes even and odd
% Gabor filters by multiplying a Gaussian with a cosine and a sine,
% respectively. It then normalizes it using the sum of squares of evenFilt
% and oddFilt.
%%
function [evenFilt, oddFilt] = gabor_filter(x, sig, sf)
    evenFilt = exp(-(x.^2)./(2*sig^2)) .* cos(2*pi*sf*x);
    oddFilt = exp(-(x.^2)./(2*sig^2)) .* sin(2*pi*sf*x);
    integral = sum(evenFilt.^2 + oddFilt.^2);
    evenFilt = evenFilt / integral;
    oddFilt = oddFilt / integral;
end