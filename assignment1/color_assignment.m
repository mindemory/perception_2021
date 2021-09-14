clear; clc; close all;

% Loading all the datasets
load surfaces;
load illuminants;
load cones;
load phosphors;

%% Q1
spectrum = linspace(400, 700, 31);

% a)
spectral_signal_18_ciea = macbeth(18, :)' .* cie_a';
figure;
plot(spectrum,spectral_signal_18_ciea,'r', 'LineWidth', 2);
xlabel('Wavelength (nm)');
ylabel('Reflected Energy');
title('Reflected energy off 18th surface illuminated by illuminant A')

% b)
spectral_signal_18_fluroescent = macbeth(18, :)' .* flourescent';
figure;
plot(spectrum,spectral_signal_18_fluroescent,'g', 'LineWidth', 2);
xlabel('Wavelength (nm)');
ylabel('Reflected Energy');
title('Reflected energy off 18th surface illuminated by fluorescent illuminant')

% c)
coneResponse_18th_ciea = cones * spectral_signal_18_ciea

% d)
coneResponse_18th_fluorescent = cones * spectral_signal_18_fluroescent

% e)
% The fluorescent light will appear more blueish as the cone response of S
% cone is higher in case of fluorescent light compared to that for the
% illuminant A. S cone responds to shorter wavelengths which are
% interpreted as blueish by the brain.

%% Q2
spectral_signal_1_ciea = macbeth(1, :)' .* cie_a';
coneResponse_1st_ciea = cones * spectral_signal_1_ciea;
cone_phosphors = cones * phosphors';
monitorSignals = inv(cone_phosphors) * coneResponse_1st_ciea