clear; close all; clc;

%% 1)
deltaT = 1; % ms
duration = 1000; % ms
t = 0:deltaT:duration-deltaT;
x = zeros(size(t));
x(100) = 1;
tau = 25; % ms

y1 = output(x, t, deltaT, tau);
fig1 = figure();
plot(t, y1, 'ro-');
xlabel('Time (ms)')
ylabel('Output')

%% a)
% Impulse response
%%
deltaT = 1; % ms
duration = 1000; % ms
t = 0:deltaT:duration-deltaT;
x = zeros(size(t));
x(100) = 1;
taus = [25, 50, 10]; % ms

fig2 = figure();
for ll = 1:length(taus)
    tau = taus(ll);
    y1 = output(x, t, deltaT, tau);
    
    subplot(1, 3, ll)
    plot(t, y1, 'r-');
    hold on;
    
    t_new = t-100;
    exponential = exp(-(t_new./tau));
    plot(exponential, 'b');
    
    xlabel('Time (ms)')
    ylabel('Output')
    ylim([0, max(y1)])
end

%% b)
% Step response
%%
deltaT = 1; % ms
duration = 1000; % ms
t = 0:deltaT:duration-deltaT;
x = zeros(size(t));
x(100:1000) = 1;
taus = [25, 50, 10]; % ms

fig3 = figure();
for ll = 1:length(taus)
    tau = taus(ll);
    y1 = output(x, t, deltaT, tau);
    
    subplot(1, 3, ll)
    plot(t, y1, 'ro', 'DisplayName', 'step resp');
    hold on;
    
    t_new = t-100;
    exponential = 1 - exp(-(t_new./tau));
    plot(exponential, 'b', 'LineWidth', 2, 'DisplayName', '1 - exp');
    
    xlabel('Time (ms)')
    ylabel('Output')
    ylim([0, max(y1)])
    legend('Location', 'southeast')
end

%% c)
% Sinusoidal response
%%
deltaT = 1; % ms
duration = 1000; % ms
t = 0:deltaT:duration-deltaT;
%x = sin(2*pi*t);
freqs = [2, 4, 8]; % ms
tau = 25;

fig4 = figure();
for ll = 1:length(freqs)
    f = freqs(ll);
    x = sin(2*pi*f*t);
    y1 = output(x, t, deltaT, tau);
    
    subplot(1, 3, ll)
    plot(t, y1, 'r-', 'DisplayName', 'step resp');
    hold on;
    
%     t_new = t-100;
%     exponential = 1 - exp(-(t_new./tau));
%     plot(exponential, 'b', 'LineWidth', 2, 'DisplayName', '1 - exp');
%     
    xlabel('Time (ms)')
    ylabel('Output')
    %ylim([0, max(y1)])
    legend('Location', 'southeast')
end

%% 2)
deltaT = 1; % ms
duration = 1000; % ms
t = 0:deltaT:duration-deltaT;
x = zeros(size(t));
x(100) = 1;
tau = 50; % ms
[f1, f2] = lp_filter(x, t, deltaT, tau);

fig5 = figure();
plot(t, f1, 'DisplayName', 'f_1')
hold on;
plot(t, f2, 'DisplayName', 'f_2')
legend()

%% 3)
deltaX = 1/120;
x_x = -2:deltaX:2;
x_y = -2:deltaX:2;
t = 0:deltaT:duration-deltaT;
x = zeros(length(x_x), length(x_y), length(t));
x(241, :, 100) = 1;
[f1, f2] = space_time_filter(x, t, deltaT, tau);

%%
sig = 0.1;
sf = 4;
evenFilt = exp(-(x_x.^2)./(2*sig^2)) .* cos(2*pi*sf*x_x);
oddFilt = exp(-(x_x.^2)./(2*sig^2)) .* sin(2*pi*sf*x_x);
integral = sum(evenFilt.^2 + oddFilt.^2);
evenFilt = evenFilt / integral;
oddFilt = oddFilt / integral;

%%
oddFast = zeros(length(x_x), length(x_y), length(t));
oddSlow = zeros(length(x_x), length(x_y), length(t));
evenSlow = zeros(length(x_x), length(x_y), length(t));
evenFast = zeros(length(x_x), length(x_y), length(t));

for tt = 1:length(t)
    oddFast(:, :, tt) = conv2(f1(:, :, tt), oddFilt', 'same');
    oddSlow(:, :, tt) = conv2(f2(:, :, tt), oddFilt', 'same');
    evenSlow(:, :, tt) = conv2(f2(:, :, tt), evenFilt', 'same');
    evenFast(:, :, tt) = conv2(f1(:, :, tt), evenFilt', 'same');
end
%%
figure();
subplot(2,2,1)
imagesc(squeeze(oddFast(:, 241, :))')
colormap(gray)

subplot(2,2,2)
imagesc(squeeze(oddSlow(:, 241, :))')
colormap(gray)

subplot(2,2,3)
imagesc(squeeze(evenFast(:, 241, :))')
colormap(gray)

subplot(2,2,4)
imagesc(squeeze(evenSlow(:, 241, :))')
colormap(gray)

%% b)
leftEven = oddFast + evenSlow;
leftOdd = -oddSlow + evenFast;
rightEven = -oddFast + evenSlow;
rightOdd = oddSlow + evenFast;

%%
figure();
subplot(2,2,1)
imagesc(squeeze(leftEven(:, 241, :))')
colormap(gray)

subplot(2,2,2)
imagesc(squeeze(leftOdd(:, 241, :))')
colormap(gray)

subplot(2,2,3)
imagesc(squeeze(rightEven(:, 241, :))')
colormap(gray)

subplot(2,2,4)
imagesc(squeeze(rightOdd(:, 241, :))')
colormap(gray)

%% c)
leftEnergy = leftEven.^2 + leftOdd.^2; 
rightEnergy = rightEven.^2. + rightOdd.^2;

figure();
subplot(2,2,1)
imagesc(squeeze(leftEnergy(:, 241, :))')
colormap(gray)
title('Left Energy')

subplot(2,2,2)
imagesc(squeeze(rightEnergy(:, 241, :))')
colormap(gray)
title('Right Energy')

%% d)



%% Functions
function y1 = output(x, t, deltaT, tau)
    y1 = zeros(length(t), 1);
    for tt = 1:length(t) - 1
        deltaY1 = (deltaT/tau) * (-y1(tt) + x(tt));
        y1(tt + 1) = y1(tt) + deltaY1;
        %y1Save(tt) = y1;
    end
end

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
    f1 = y(:, 3) - y(:, 5);
    f2 = y(:, 5) - y(:, 7);
end


function [f1, f2] = space_time_filter(x, t, deltaT, tau)
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
        f1(:, :, tt) = y(:, :, tt, 3) - y(:, :, tt, 5);
        f2(:, :, tt) = y(:, :, tt, 5) - y(:, :, tt, 7);
    end
    
end
% function space_time_filter(x, t, deltaT, tau)
%     n = [3,5,5,7];
%     for tt = 1:size(x,1)
%         % Temporal filters
%         deltaY = (deltaT/tau) * (- y(1,:) + input(tt,:));
%         y(1,:) = y(1,:) + deltaY;
%         for nn = 2:max(n)
%             deltaY = (deltaT/tau) * (-y(nn,:) + y(nn-1,:));
%             y(nn,:) = y(nn,:) + deltaY; end
%         rtFast = y(n(1),:)-y(n(2),:);
%         rtSlow = y(n(3),:)-y(n(4),:);
%         % Spatial filters
%         oddFast = spatialConvolution(rtFast,oddFilt);
%         oddSlow = spatialConvolution(rtSlow,oddFilt);
%         evenSlow = spatialConvolution(rtSlow,evenFilt);
%         evenFast = spatialConvolution(rtFast,evenFilt);
%         % Direction selective filters and motion energy
%         leftEven = oddFast + evenSlow;
%         leftOdd = -oddSlow + evenFast;
%         leftEnergy = leftEven.^2 + leftOdd.^2; rightEven = -oddFast + evenSlow;
%         rightOdd = oddSlow + evenFast;
%         rightEnergy = rightEven.^2. + rightOdd.^2;
%     end
% end