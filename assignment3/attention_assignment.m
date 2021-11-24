clear; close all; clc;

%% 1)
% The role of attention is to increase the discriminability of the stimulus
% and suppress the unwanted noise. Therefore, in terms of spatial
% frequency, I would expect a broadening of spatial frequency and
% orientation bandwidths because of attention. This would ensure that the
% sensitivity of all the spatial frequencies and all orientation bandwidths
% is increased due to attention. As a result, the performance of the
% subject would improve for all spatial frequencies and orientation
% bandwidths. However, this would also imply an increased sensitivity to
% noise at all spatial frequencies and orientation bandwidths.
%%
% According to this hypothesis, the effect of attention can be modeled as
% an increase in the Gaussian response curve over all spatial frequencies
% and orientation bandwidths. This can be modeled using Naka-Rushton
% functions and their modifications as in (Wang et al., 2017; Jigo and
% Carrasco, 2020). However, for the sake of this assigment, I have kept the
% approximation to a Gaussian model that serves the purpose. In fact the
% cdf of a Gaussian model is a psychometric curve that can model both
% response gain and/or frequency and/or orientation gains.
%%
% In terms of a psychological experiment, attention can be modeled as a
% cue. And absense or no attention can be modeled by a neutral cue. At the
% same time, the location that is uncued will serve as unattended stimulus.
% Modeling this using a Gaussian for spatial frequency and orientation of
% the stimulus, would provide following curves in the case that the
% hypothesis does hold true.
%%
mu_sf = 3;
sigma_sf = 1;
sf = mu_sf - 0.5 * sigma_sf:0.01:mu_sf+1.5*sigma_sf;
sens_sf_neutral = normpdf(sf, mu_sf, sigma_sf);
sens_sf_cued = normpdf(sf, mu_sf, sigma_sf) * 1.1;
sens_sf_uncued = normpdf(sf, mu_sf, sigma_sf) / 1.1;

ori = -90 : 1: 90;
mu_ori = 0;
sigma_ori = 40;
sens_ori_neutral = normpdf(ori, mu_ori, sigma_ori);
sens_ori_cued = normpdf(ori, mu_ori, sigma_ori) * 1.2;
sens_ori_uncued = normpdf(ori, mu_ori, sigma_ori) / 1.2;

fig1 = figure();
suptitle(sprintf('Effect of attention on spatial frequency and orientation\n'))

subplot(1, 2, 1)
plot(sf, sens_sf_neutral, 'DisplayName', 'neutral')
hold on;
plot(sf, sens_sf_cued, 'DisplayName', 'cued')
plot(sf, sens_sf_uncued, 'DisplayName', 'uncued')
xlabel('Spatial frequency (cpd)')
ylabel('Sensitivity (AU)')
legend();

subplot(1, 2, 2)
plot(ori, sens_ori_neutral, 'DisplayName', 'neutral')
hold on;
plot(ori, sens_ori_cued, 'DisplayName', 'cued')
plot(ori, sens_ori_uncued, 'DisplayName', 'uncued')
xlabel('Orientation (degrees)')
ylabel('Sensitivity (AU)')
xlim([-90, 90])
legend();

%%
% However, the role of attention would also depend on whether endogenous or
% exogenous attention are deployed, the demands of the task, and the nature 
% of the task (Carrasco and Yeshurun, 2009; Barbot and Carrasco, 2017;
% Jigo and Carrasco, 2020; Dugue et al., 2020). In these and other
% instances, there are several similarities in behavioral performance for
% endogenous covert and exogenous covert attention despite quite a few
% dissimilarities based on task demand. Here, I am assuming that the
% performance for both endogenous and exogenous would be same or close to
% similar. An experiment can be run to verify if there are differences in
% the effect of endogenous vs exogenous attention on spatial frequency and
% orientation bandwidths.
%%
% Adapting the task design from Jigo and Carrsco, 2020, I have designed a
% behavioral experiment that would help determine the sensitivity of the
% subject to different spatial frequencies and orientations based on
% whether the stimulus is attended (cued), not attended (neutral) or
% unattended (uncued). Additionaly, the experimental setup can also help
% determine the differences (if any) between the performance in an
% exogenous vs endogenous attention task.
%%
%% Exogenous attention task



% b1 = normcdf(sf, mu_sf, sigma_sf);
% b2 = normcdf(sf, mu_ori, sigma_ori) * 2;
% fig2 = figure();
% plot(sf, b1, 'DisplayName', 'a1')
% hold on;
% plot(sf, b2, 'DisplayName', 'a2')
% legend();

% function df = naka_rushton(d_max, c, n, f, condition)
%     if condition == 'neutral'
%         c50_f = c50N(f);
%     elseif condition == 'cued'
%         c50_f = c50A(f);
%     end
%     df = d_max * (c^n / (c^n + c50_f^n));
%     
% end
% 
% function c50_f = c50N(f)
%     m = 2;
%     f_csf = 1;
%     s = 0.5;
%     c50_f = 1/(m * f^(f_csf/s) * e^(-f/s));
% end
% 
% function c50_f = c50A(f)
%     p = 2;
%     sigma = 0.5;
%     exponent = (f - f_benefit)/sig;
%     b = gamma * e^(-exponent^p) + delta;
%     c50_f = c50N(f)/b;
% end