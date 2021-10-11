clear; close all; clc;

%% Q1a)
stim_dur = linspace(0.1, 2, 20);
trials_per_experiment = 100;
contrast = 0.2;
var = 1;

trialspace = rand(trials_per_experiment, 1);
trialspace = (trialspace>.5) -1 * (trialspace<=.5);
percentage_correct_simul = [];
dprime_simul = [];
percentage_correct_theor = [];
dprime_theor = [];
for dur = stim_dur
    time_steps = dur * 10;

    B_mean = -2 * time_steps * contrast;
    A_mean = 2 * time_steps * contrast;
    B_var = time_steps * var;
    A_var = time_steps * var;
    optim_criterion = (A_mean + B_mean) / 2;
    
    evidence = time_steps * 2 * contrast * trialspace + ...
        trialspace .* randn(trials_per_experiment, 1);

    response = (evidence > optim_criterion) -1 * (evidence <= optim_criterion);
    H = 0;
    M = 0;
    CR = 0;
    FA = 0;
    for i = 1:trials_per_experiment
        if trialspace(i) == 1
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
    
    pcorrect = (H + CR) * 100/trials_per_experiment;
    percentage_correct_simul = [percentage_correct_simul pcorrect];
    Hrate = H/(H+M);
    FArate = FA/(FA+CR);
    zH = -sqrt(2) * erfcinv(2*Hrate);
    zFA = -sqrt(2) * erfcinv(2*FArate);
    dp = zH - zFA;
    dprime_simul = [dprime_simul dp];
    
    
    

    resp_range = linspace(B_mean-2*B_var, A_mean+2*A_var, 100);
    B_pdf = normpdf(resp_range, B_mean, sqrt(B_var));
    B_cdf = normcdf(resp_range, B_mean, sqrt(B_var));
    A_pdf = normpdf(resp_range, A_mean, sqrt(A_var));
%     figure();
%     plot(resp_range, [B_pdf', A_pdf'])
%     ylabel('Probability')
%     xlabel('Total evidence')
%     legend('B', 'A')
%     title('A and B pdfs')
    percentage_correct_theor = [percentage_correct_theor, sum(A_pdf .* B_cdf)];
    dprime_theor = [dprime_theor, A_mean - B_mean];
end
figure();
plot(stim_dur, percentage_correct_simul, 'b-o', 'LineWidth', 2, ...
    'MarkerSize', 5)
xlabel('Duration of stimulus (s)')
ylabel('Percentage correct')
figure();
plot(stim_dur, dprime_simul, 'm-o', 'LineWidth', 2, ...
    'MarkerSize', 5);
xlabel('Duration of stimulus (ms)')
ylabel("d'")

figure();
plot(stim_dur, percentage_correct_theor, 'b-o', 'LineWidth', 2, ...
    'MarkerSize', 5)
xlabel('Duration of stimulus (s)')
ylabel('Percentage correct')
figure();
plot(stim_dur, dprime_theor, 'm-o', 'LineWidth', 2, ...
    'MarkerSize', 5);
xlabel('Duration of stimulus (ms)')
ylabel("d'")
