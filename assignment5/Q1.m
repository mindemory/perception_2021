clear; close all; clc;

%% a)
% There are two depth cues for the distance to the beach. Both of the cues
% are independent Gaussian random variables. These are unbiased cues
% implies that the expected value of the cue is equal to the actual value
% of the depth.

%% b)
% Let the depth estimate due to the first cue be $$D_1 $$ and the depth
% estimate due to the second cue be $$D_2 $$. Let the weightage given to
% the first cue be w. Therefore, the weightage given to the second cue will
% be 1-w. Hence, the combined estimate of depth becomes:
%%
% $$D = wD_1 + (1 - w)D_2 $$
%%
% Therefore, the variance of the estimate will be:
%%
% $$Var(D) = Var(wD_1 + (1 - w)D_2) = Var(wD_1)+ Var((1 - w)D_2) $$
%%
% Therefore,
%%
% $$Var(D) = w^2 Var(D_1)+ (1 - w)^2 Var(D_2) $$
%%
% An optimal observer would try to minimize the resulting estimate of the
% variance of depth. Therefore, the goal of the optimal observer would be
% to make the derivative of Var(D) = 0.
%%
% $$\frac{d Var(D)}{dw} = 2w Var(D_1)- 2(1 - w) Var(D_2) = 0$$
%%
% Solving for w:
%%
% $$w Var(D_1)- (1 - w) Var(D_2) = 0 $$
%%
% Therefore,
%%
% $$w Var(D_1)- Var(D_2) + w Var(D_2) = 0 $$
%%
% Therefore,
%%
% $$w (Var(D_1) + Var(D_2)) = Var(D_2) $$
%%
% Therefore,
%%
% $$w = \frac{Var(D_2)}{Var(D_1) + Var(D_2)} $$
%%
% We can check by taking the second-derivative of V that this is indeed the
% w that minimizes the variance of the estimate.
%%
% The given variances of the two cues are:
%%
% $$Var(D_1) = 4, Var(D_2) = 1 $$
%%
% Substituting back into the formula for w, we get:
%%
% $$w = \frac{1}{4 + 1} = \frac{1}{5} = 0.2 $$
%%
% Also,
%%
% $$1 - w = 1 - \frac{1}{5} = \frac{4}{5} = 0.8 $$
%%
% Therefore, an optimal observer would give a weight of 0.2 to the first
% cue and a weight of 0.8 to the second cue. Plugging the weights back into
% the formula for computing the variance of the output, we get:
%%
% $$Var(D) = 0.2^2\times 4+ 0.8^2\times 1 $$
%%
% Therefore,
%%
% $$Var(D) = 0.04\times 4+ 0.64\times 1 $$
%%
% Therefore,
%%
% $$Var(D) = 0.16+ 0.64 $$
%%
% Therefore,
%%
% $$Var(D) = 0.80 $$

%% c)
% If, instead, equal weigths were given to the two cues, then we have:
%%
% $$w = 1 - w = 0.5 $$
%%
% Therefore, the estimates of the variance of the depth now becomes:
%%
% $$Var(D) = 0.5^2\times 4+ 0.5^2\times 1 $$
%%
% Therefore,
%%
% $$Var(D) = 0.25\times 4+ 0.25\times 1 $$
%%
% Therefore,
%%
% $$Var(D) = 1+ 0.25 $$
%%
% Therefore,
%%
% $$Var(D) = 1.25 $$
%%
% The difference between the variance of the optimal estimator and this estimator is
% 1.25 - 0.8 = 0.45
%% d)
% As can be seen, the "average cue" has a higher variance than the second
% cue. Therefore, David would be better off selecting just the second cue
% and ignoring the first cue instead of averaging the two cues.