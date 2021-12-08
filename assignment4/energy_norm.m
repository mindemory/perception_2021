%% energy normalization
% The function computes normalized energies for each neuron by normalizing
% it with the sum of energies of responses for all the neurons. stdev here
% controls the normalization.
%%
function [leftEnergyNorm, rightEnergyNorm, upEnergyNorm, downEnergyNorm] = ...
    energy_norm(leftEnergy, rightEnergy, upEnergy, downEnergy)
    sumEnergy = leftEnergy + rightEnergy + upEnergy + downEnergy;
    stdev = 0.02;
    leftEnergyNorm = leftEnergy ./ (sumEnergy + stdev^2);
    rightEnergyNorm = rightEnergy ./ (sumEnergy + stdev^2);
    upEnergyNorm = upEnergy ./ (sumEnergy + stdev^2);
    downEnergyNorm = downEnergy ./ (sumEnergy + stdev^2);
end
