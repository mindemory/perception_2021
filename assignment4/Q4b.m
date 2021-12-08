clear; close all; clc;

%% b)
% Repeating the experiment for a cross-orientation experiment where the
% contrast of the rightward moving grating is varied between 1% and 50% and
% the contrast of the upward moving grating is kept at 50%. The normalized
% energies of each of the direction-selective neurons are computed for
% different contrasts of the rightward moving gratins when the upward
% moving grating is either present or absent. The rightward-selective
% neuron responds the strongest when there is only rightward moving grating
% and the response is dampened by the presence of the upward selective
% grating. However, the effect of upward grating decreases with an increase
% in the contrast of the rightward moving grating. A similar pattern is
% observed for the leftward selective response. However, the scale of the
% responses for the leftward selective neuron are lower than the rightward
% selective neuron as the stimulus is not along the preferred direction of
% the neuron. For the upward selective neuron, the response is higher when
% the contrast of the rightward moving grating is lower and the responses
% decrease when the contrast of the rightward moving grating increases.
% Similar pattern is observed for the downward selective neuron but the
% responses are lower than the upward selective neuron as the stimulus is
% not along the preferred direction of the neuron.
%%
tic
deltaX = 1/120; % spatial sampling rate
x_x = -2:deltaX:2; % spatial array along x axis
x_y = -2:deltaX:2; % spatial array along y axis
deltaT = 1; % ms
duration = 1000; % ms

t = 0:deltaT:duration-deltaT; % time-array

tau = 25; % ms

phase_shift = 2*pi/125;
contrasts = [1, 5, 10, 20, 35, 50];
cross_oris = [0, 1];
phase = 0; % initial phase of the stimulus
sf = 30; % cycles/pixel
sig = 0.1; % standard deviation of the Gaussian (in deg)

[evenFilt, oddFilt] = gabor_filter(x_x, sig, 4); % computes spatial filters

ori = "lr"; % for left-right

mean_energies = zeros(length(contrasts), 4, length(cross_oris));

for ff = 1:2
    cross_ori = cross_oris(ff);
    for ii = 1:length(contrasts)
        contrast = contrasts(ii);
        [leftEven, leftOdd, rightEven, rightOdd, leftEnergy, rightEnergy, ...
                upEven, upOdd, downEven, downOdd, upEnergy, downEnergy] = ...
                neuron_responses_cross_ori(x_x, x_y, t, deltaT, tau, contrast, phase, ...
                phase_shift, sf, oddFilt, evenFilt, cross_ori);
        [leftEnergyNorm, rightEnergyNorm, upEnergyNorm, downEnergyNorm] = ...
            energy_norm(leftEnergy, rightEnergy, upEnergy, downEnergy);
        x_y_dim = 241;
        st_tm = 400; % ms
        mean_energies(ii, :, ff) = [mean(leftEnergyNorm(x_y_dim, x_y_dim, st_tm:end)), ...
            mean(rightEnergyNorm(x_y_dim, x_y_dim, st_tm:end)), ...
            mean(upEnergyNorm(x_y_dim, x_y_dim, st_tm:end)), ...
            mean(downEnergyNorm(x_y_dim, x_y_dim, st_tm:end))];
    end
end
toc

plt_titles = ["leftEnergyNorm mean", "rightEnergyNorm mean", ...
    "upEnergyNorm mean", "downEnergyNorm mean"];

figure()
for ss = 1:size(mean_energies, 2)
    subplot(2, 2, ss)
    plot(contrasts, mean_energies(:, ss, 1), 'o-', ...
        'DisplayName', 'right', 'LineWidth', 2);
    legend('location', 'northwest')
    hold on;
    plot(contrasts, mean_energies(:, ss, 2), 'o-', ...
        'DisplayName', 'up+right', 'LineWidth', 2);
    xlabel('Contrast of rightward grating')
    ylabel('Normalized Energy')
    title(plt_titles(ss))
end
