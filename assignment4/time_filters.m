%%  time filters
% Similar to lp_filter but instead computes the temporal filters over the
% space.
%%
function [f1, f2] = time_filters(x, t, deltaT, tau)
    
    [x_size, y_size, t_length] = size(x);
    y = zeros(x_size, y_size, t_length, 7);
    f1 = zeros(x_size, y_size, t_length);
    f2 = zeros(x_size, y_size, t_length);

    for tt = 1:t_length - 1
        for jj = 1:7
            if jj == 1
                deltaY = (deltaT / tau) * (-y(:, :, tt, jj) + x(:, :, tt));
                y(:, :, tt + 1, jj) = y(:, :, tt, jj) + deltaY;
            else
                deltaY = (deltaT / tau) * (-y(:, :, tt, jj) + y(:, :, tt, jj - 1));
                y(:, :, tt + 1, jj) = y(:, :, tt, jj) + deltaY;
            end
        end
        f1(:, :, tt) = y(:, :, tt, 3) - y(:, :, tt, 5); % Fast filter
        f2(:, :, tt) = y(:, :, tt, 5) - y(:, :, tt, 7); % Slow filter
    end
end
