clear; close all; clc;

%%
% lp_filter uses a cascade of exponential low-pass filters and computes the
% temporal filters f1 and f2
%%
deltaT = 1; % ms
duration = 1000; % ms
t = 0:deltaT:duration-deltaT;
x = zeros(size(t));
x(100) = 1;
tau = 25; % ms
[f1, f2] = lp_filter(x, t, deltaT, tau);

figure();
plot(t, f1, 'DisplayName', 'f_1')
hold on;
plot(t, f2, 'DisplayName', 'f_2')
xlabel('time (ms)')
title('Temporal filters')
legend()
