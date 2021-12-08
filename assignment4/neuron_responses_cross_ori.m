%% neuron response for cross-orientation
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
function [leftEven, leftOdd, rightEven, rightOdd, leftEnergy, rightEnergy, ...
            upEven, upOdd, downEven, downOdd, upEnergy, downEnergy] = ...
            neuron_responses_cross_ori(x_x, x_y, t, deltaT, tau, contrast, phase, ...
            phase_shift, sf, oddFilt, evenFilt, cross_ori)
    
    stim = zeros(length(x_x), length(x_y), length(t)); % Initializing the stimulus
    
    amp_sinusoid_right = contrast * 1/100; % amplitude of the rightward moving sinusoid is scaled by the contrast of the stimulus desired
    amp_sinusoid_up = 0.5; % amplitude of the upward moving sinusoid is 50%
    phase_up = phase; % Initializing phase of upward moving sinusoid
    phase_right = phase; % Initializing phase of rightward moving sinusoid
    if cross_ori == 1
    
        for tt = 1:length(t)
            phase_right = phase_right - phase_shift; % Shifting the phase at each time-step to move the gratings rightward
            phase_up = phase_up + phase_shift; % Shifting the phase at each time-step to move the gratings upward
            rightwardsinusoid = get_grating(x_x, amp_sinusoid_right, phase_right,...
                sf, "lr"); % Generating rightward moving 2-D sinusoid at each time-step
            upwardsinusoid = get_grating(x_x, amp_sinusoid_up, phase_up,...
                sf, "ud"); % Generating upward moving 2-D sinusoid at each time-step
            stim(:, :, tt) = rightwardsinusoid + upwardsinusoid; % Stimulus is the sum of the two gratings
        end
        
    else
        for tt = 1:length(t)
            phase_right = phase_right - phase_shift; % Shifting the phase at each time-step to move the gratings rightward
            rightwardsinusoid = get_grating(x_x, amp_sinusoid_right, phase_right,...
                sf, "lr"); % Generating rightward moving 2-D sinusoid at each time-step
            stim(:, :, tt) = rightwardsinusoid; % Stimulus is the rightward moving grating
        end

    end
    [f1, f2] = time_filters(stim, t, deltaT, tau); % Obtain temporal filters with sinusoid as input
        
    [oddFastlr, oddSlowlr, evenFastlr, evenSlowlr] = ...
        conv_filts(f1, f2, oddFilt, evenFilt); % convolve temporal and horizontal spatial filters
    [oddFastud, oddSlowud, evenFastud, evenSlowud] = ...
        conv_filts(f1, f2, oddFilt', evenFilt'); % convolve temporal and vertical spatial filters
    
    [leftEven, leftOdd, rightEven, rightOdd, leftEnergy, rightEnergy] = ...
        motion_energy(oddFastlr, oddSlowlr, evenFastlr, evenSlowlr); % compute even, odd, energy for leftward and rightward oriented neurons
    [upEven, upOdd, downEven, downOdd, upEnergy, downEnergy] = ...
        motion_energy(oddFastud, oddSlowud, evenFastud, evenSlowud); % compute even, odd, energy for upward and downward oriented neurons
    
end        