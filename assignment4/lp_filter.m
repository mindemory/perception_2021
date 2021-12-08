%%  low-pass filter
% Given the input x, time_array t, time-step deltaT, and time-constant
% tau, the function solves a cascade of ODEs to compute the temporal
% filters f1 and f2. Here f1 is a fast-filter while f2 is a slow filter.
%%
function [f1, f2] = lp_filter(x, t, deltaT, tau)
    
    y = zeros(length(t), 7);
    
    for tt = 1:length(t) - 1
        for jj = 1:7
            if jj == 1
                deltaY = (deltaT / tau) * (-y(tt, jj) + x(tt));
                y(tt + 1, jj) = y(tt, jj) + deltaY;
            else
                deltaY = (deltaT / tau) * (-y(tt, jj) + y(tt, jj - 1));
                y(tt + 1, jj) = y(tt, jj) + deltaY;
            end
        end
    end
    f1 = y(:, 3) - y(:, 5); % Fast filter
    f2 = y(:, 5) - y(:, 7); % Slow filter
end