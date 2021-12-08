%% plotting 
% The function creates 4-subplots: 1 for each input passed. Based on the
% orientation of the temporal filter used in computing the inputs, it will
% create a x-t/y-t slice for visualization. If, however, the inputs are
% energies, it will plot the left and right energies by taking x-t
% slices and the up and down energies by taking y-t slices.
%%
function plot_func(a, b, c, d, a_title, b_title, c_title, d_title, direc)
    
    if strcmp(direc, 'lr') % for left-right
        a = squeeze(a(241, :, :))';
        b = squeeze(b(241, :, :))';
        c = squeeze(c(241, :, :))';
        d = squeeze(d(241, :, :))';
    elseif strcmp(direc, 'ud') % for up-down
        a = squeeze(a(:, 241, :))';
        b = squeeze(b(:, 241, :))';
        c = squeeze(c(:, 241, :))';
        d = squeeze(d(:, 241, :))';
    elseif strcmp(direc, 'energy')
        a = squeeze(a(241, :, :))';
        b = squeeze(b(241, :, :))';
        c = squeeze(c(:, 241, :))';
        d = squeeze(d(:, 241, :))';
    end

    subplot(2,2,1)
    imagesc(a)
    colormap(gray)
    xticks([1, 121, 241, 361, 481])
    xticklabels([-2, -1, 0, 1, 2])
    yticks(0:100:1000)
    xlabel('Visual angle (deg)')
    ylabel('time (ms)')
    title(a_title)

    subplot(2,2,2)
    imagesc(b)
    colormap(gray)
    xticks([1, 121, 241, 361, 481])
    xticklabels([-2, -1, 0, 1, 2])
    yticks(0:100:1000)
    xlabel('Visual angle (deg)')
    ylabel('time (ms)')
    title(b_title)

    subplot(2,2,3)
    imagesc(c)
    colormap(gray)
    xticks([1, 121, 241, 361, 481])
    xticklabels([-2, -1, 0, 1, 2])
    yticks(0:100:1000)
    xlabel('Visual angle (deg)')
    ylabel('time (ms)')
    title(c_title)

    subplot(2,2,4)
    imagesc(d)
    colormap(gray)
    xticks([1, 121, 241, 361, 481])
    xticklabels([-2, -1, 0, 1, 2])
    yticks(0:100:1000)
    xlabel('Visual angle (deg)')
    ylabel('time (ms)')
    title(d_title)
end
