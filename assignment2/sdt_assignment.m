clear; close all; clc;

%% Q1a)
%% variable duration and constant contrast
stim_dur = 0.1: 0.1: 2;
trials_per_experiment = 10000;
contrast = 0.2;
var = 1;

trialspace = rand(trials_per_experiment, 1);
trialspace = (trialspace>.5) -1 * (trialspace<=.5);
percentage_correct_simul = zeros(size(stim_dur));
dprime_simul = zeros(size(stim_dur));
percentage_correct_theor = zeros(size(stim_dur));
dprime_theor = zeros(size(stim_dur));
crit = 'None';
for i = 1:length(stim_dur)
    dur = stim_dur(i);
    [optim_criterion, pcorr_simul, dp_simul, pcorr_theor, dp_theor] = ...
        run_experiment(dur, contrast, var,...
        trialspace, trials_per_experiment, crit);
    if i == 1
        disp(['The optimum criterion is: ', num2str(optim_criterion)])
    end
    percentage_correct_simul(i) = pcorr_simul;
    dprime_simul(i) = dp_simul;
    percentage_correct_theor(i) = pcorr_theor;
    dprime_theor(i) = dp_theor;
end
figure(1);
xlab = 'Duration of stimulus (s)';
plot_figs(stim_dur, percentage_correct_simul, dprime_simul, percentage_correct_theor, ...
    dprime_theor, xlab)

%% constant duration and variable contrast
dur = 0.2;
stim_contrast = 0.05: 0.05: 1;

percentage_correct_simul = zeros(size(stim_contrast));
dprime_simul = zeros(size(stim_contrast));
percentage_correct_theor = zeros(size(stim_contrast));
dprime_theor = zeros(size(stim_contrast));
for i = 1:length(stim_contrast)
    contrast = stim_contrast(i);
    [optim_criterion, pcorr_simul, dp_simul, pcorr_theor, dp_theor] = ...
        run_experiment(dur, contrast, var,...
        trialspace, trials_per_experiment, crit);
    if i == 1
        disp(['The optimum criterion is: ', num2str(optim_criterion)])
    end
    percentage_correct_simul(i) = pcorr_simul;
    dprime_simul(i) = dp_simul;
    percentage_correct_theor(i) = pcorr_theor;
    dprime_theor(i) = dp_theor;
end

figure(2);
xlab = 'Contrast of stimulus (s)';
plot_figs(stim_contrast, percentage_correct_simul, dprime_simul, percentage_correct_theor, ...
    dprime_theor, xlab)


%%%%%%%%
%% Q1b)
dur = 0.4;
trials_per_experiment = 100;
contrast = 0.1;
var = 1;

trialspace = rand(trials_per_experiment, 1);
trialspace = (trialspace>.25) -1 * (trialspace<=.25);

[optim_criterion, pcorr_simul, dp_simul, pcorr_theor, dp_theor] = ...
    run_experiment(dur, contrast, var,...
    trialspace, trials_per_experiment, crit);
disp(['The optimum criterion is: ', num2str(optim_criterion)])

%figure(3);
%xlab = 'Contrast of stimulus (s)';
%plot_figs(stim_contrast, percentage_correct_simul, dprime_simul, percentage_correct_theor, ...
%    dprime_theor, xlab)


function [optim_criterion, pcorr_simul, dp_simul, pcorr_theor, dp_theor] = ...
    run_experiment(dur, contrast, var, trials, t_length, crit)
    %% Simulation
    time_steps = dur * 10;
    B_mean = -2 * time_steps * contrast;
    A_mean = 2 * time_steps * contrast;
    B_var = time_steps * var;
    A_var = time_steps * var;
    optim_criterion = (A_mean + B_mean) / 2;
    if crit == 'None'
        crit = optim_criterion;
    end
    
    evidence = time_steps * 2 * contrast * trials + ...
        sqrt(A_var) * trials .* randn(t_length, 1);

    response = (evidence > crit) -1 * (evidence <= crit);
    H = 0; M = 0; CR = 0; FA = 0;
    for i = 1:t_length
        if trials(i) == 1
            if response(i) == 1
                H = H + 1;
            else
                M = M + 1;
            end
        else
            if response(i) == -1
                CR = CR + 1;
            else
                FA = FA + 1;
            end
        end
    end
    
    pcorr_simul = (H + CR) * 100/t_length;
    
    Hrate = H/(H+M);
    FArate = FA/(FA+CR);
    if Hrate == 1
        Hrate = 0.9999;
    end
    if FArate == 0
        FArate = 0.0001;
    end
    zH = -sqrt(2) * erfcinv(2*Hrate);
    zFA = -sqrt(2) * erfcinv(2*FArate);
    dp_simul = zH - zFA;
    
    %% Theoretical
    resp_range = B_mean-3*B_var: 0.01: A_mean+3*A_var;
    B_cdf = normcdf(resp_range, B_mean, sqrt(B_var));
    A_pdf = normpdf(resp_range, A_mean, sqrt(A_var));
    pcorr_theor = sum(A_pdf .* B_cdf);
    dp_theor = (A_mean - B_mean)/sqrt(A_var);
end

function plot_figs(x, y1, y2, y3, y4, xlab)
    subplot(2, 2, 1)
    plot(x, y1, 'b-o', 'LineWidth', 2, ...
        'MarkerSize', 5)
    xlabel(xlab)
    ylabel('Percentage correct')
    title('Psychometric function simulation')
    subplot(2, 2,2)
    plot(x, y2, 'm-o', 'LineWidth', 2, ...
        'MarkerSize', 5);
    xlabel(xlab)
    ylabel("d'")
    title("d' simulation")

    subplot(2, 2, 3)
    plot(x, y3, 'b-o', 'LineWidth', 2, ...
        'MarkerSize', 5)
    xlabel(xlab)
    ylabel('Percentage correct')
    title('Psychometric function theoretical')
    subplot(2, 2, 4)
    plot(x, y4, 'm-o', 'LineWidth', 2, ...
        'MarkerSize', 5);
    xlabel(xlab)
    ylabel("d'")
    title("d' theoretical")
end