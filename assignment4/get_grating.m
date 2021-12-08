%% create 2-D sinusoid
% The function creates a 2-D sinusoid with given amplitude, phase, spatial
% frequency (sf), orientation ('lr' for horizontal, and 'ud' for vertical),
% and the range of x and y for the sinusoid is determined by a
% one-dimensional array x_x.
%%
function sinewave2D = get_grating(x_x, amp, phase, sf, ori)
    [X,Y] = meshgrid(x_x .* sf);
    if strcmp(ori, 'lr')
        sinewave2D = amp * sin(X + phase);
    elseif strcmp(ori, 'ud')
        sinewave2D = amp * sin(Y + phase);
    end
end
