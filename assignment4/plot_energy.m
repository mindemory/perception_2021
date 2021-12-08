%% plotting energies
% The function creates 4-subplots, each corresponding to the preferred
% direction of the neuron. Each subplot has 3 lines: even, odd and
% energy. For each line-plot, a time-slice is taken for the center of the
% screen with (x, y) = (0 deg, 0 deg)
%%
function plot_energy(t, leftEnergy, leftEven, leftOdd, rightEnergy, rightEven, ...
    rightOdd, upEnergy, upEven, upOdd, downEnergy, downEven, downOdd, title_st)
    
    sgtitle(title_st)
    pp = 241; % the x and y co-ords for 0 degree spatial angle
    
    subplot(2, 2, 1)
    plot(t, squeeze(leftEnergy(pp, pp, :)), 'DisplayName', 'leftEnergy');
    hold on;
    plot(t, squeeze(leftEven(pp, pp, :)), 'DisplayName', 'leftEven');
    plot(t, squeeze(leftOdd(pp, pp, :)), 'DisplayName', 'leftOdd');
    xlabel('Time (ms)')
    ylabel('Response')
    legend()
    
    subplot(2, 2, 2)
    plot(t, squeeze(rightEnergy(pp, pp, :)), 'DisplayName', 'rightEnergy');
    hold on;
    plot(t, squeeze(rightEven(pp, pp, :)), 'DisplayName', 'rightEven');
    plot(t, squeeze(rightOdd(pp, pp, :)), 'DisplayName', 'rightOdd');
    xlabel('Time (ms)')
    ylabel('Response')
    legend()
    
    subplot(2, 2, 3)
    plot(t, squeeze(upEnergy(pp, pp, :)), 'DisplayName', 'upEnergy');
    hold on;
    plot(t, squeeze(upEven(pp, pp, :)), 'DisplayName', 'upEven');
    plot(t, squeeze(upOdd(pp, pp, :)), 'DisplayName', 'upOdd');
    xlabel('Time (ms)')
    ylabel('Response')
    legend()
    
    subplot(2, 2, 4)
    plot(t, squeeze(downEnergy(pp, pp, :)), 'DisplayName', 'downEnergy');
    hold on;
    plot(t, squeeze(downEven(pp, pp, :)), 'DisplayName', 'downEven');
    plot(t, squeeze(downOdd(pp, pp, :)), 'DisplayName', 'downOdd');
    xlabel('Time (ms)')
    ylabel('Response')
    legend()
end