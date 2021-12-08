%%  iir_filter
% Given the input x, time_array t, time-step deltaT, and time-constant
% tau, the function solves the ODE to compute the response obtained
% throught the IIR-filter
%%
function y1 = iir_filter(x, t, deltaT, tau)
    y1 = zeros(length(t), 1);
    for tt = 1:length(t) - 1
        deltaY1 = (deltaT/tau) * (-y1(tt) + x(tt));
        y1(tt + 1) = y1(tt) + deltaY1;
    end
end

