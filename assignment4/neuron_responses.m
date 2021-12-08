%% neuron response
% The neuron_responses runs the entire experiment by computing the
% responses of leftward, rightward, upward, downward oriented neurons for a
% given series of gratings over time. The gratings are designed using
% get_grating function that creates a sinusoidal grating in space with the
% spatial frequency of 4 cyc/deg. Since the spatial scale of the space
% varies from -2 deg to 2 deg in steps of 1/120 degree. Hence, 4 degrees
% correspond to 481 pixels. Therefore, the spatial frequency of the
% sinusoid is 4 cyc/ 481 pixels. Therefore, the time-period of the sinusoid
% becomes 481/4 ~ 30 pixels/cyc. The sinusoidal grating also moves in time
% with a frequency of 8Hz = 8 cycles/sec. 1 sec = 1000 ms = 1000 frames.
% Therefore, in each frame, the sinusoid moves by 1000/8 = 125
% frames/cycle. Hence the phase of the sinusoid will change at each
% time-step by 2*pi/125 in the preferred direction. The sinusoid is then
% used as an input to compute time_filters which are convolved with the
% spatial filters to produce neuronal responses and energies.
%%
function [leftEven, leftOdd, rightEven, rightOdd, leftEnergy, rightEnergy, ...
            upEven, upOdd, downEven, downOdd, upEnergy, downEnergy] = ...
            neuron_responses(x_x, x_y, t, deltaT, tau, contrast, phase, ...
            phase_shift, sf, ori, oddFilt, evenFilt)
    
    sinusoid_x = zeros(length(x_x), length(x_y), length(t)); % Initializing the sinusoid
    
    amp_sinusoid = contrast * 1/100; % amplitude of the sinusoid is scaled by the contrast of the stimulus desired
    
    for tt = 1:1000
        phase = phase + phase_shift; % Shifting the phase at each time-step to move the gratings
        sinusoid_x(:, :, tt) = get_grating(x_x, amp_sinusoid, phase, sf, ori); % Generating 2-D sinusoid at each time-step
    end
    
    [f1, f2] = time_filters(sinusoid_x, t, deltaT, tau); % Obtain temporal filters with sinusoid as input
    
    [oddFastlr, oddSlowlr, evenFastlr, evenSlowlr] = ...
        conv_filts(f1, f2, oddFilt, evenFilt); % convolve temporal and horizontal spatial filters
    [oddFastud, oddSlowud, evenFastud, evenSlowud] = ...
        conv_filts(f1, f2, oddFilt', evenFilt'); % convolve temporal and vertical spatial filters
    
    [leftEven, leftOdd, rightEven, rightOdd, leftEnergy, rightEnergy] = ...
        motion_energy(oddFastlr, oddSlowlr, evenFastlr, evenSlowlr); % compute even, odd, energy for leftward and rightward oriented neurons
    [upEven, upOdd, downEven, downOdd, upEnergy, downEnergy] = ...
        motion_energy(oddFastud, oddSlowud, evenFastud, evenSlowud); % compute even, odd, energy for upward and downward oriented neurons
end        
