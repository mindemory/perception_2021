function renderSurfaces(surfaces,light,cones,phosphors,range,header)
%
% renderSurfaces(surfaces,light,cones,phosphors,[range],[header])
%
% Displays a matrix of tiles based on inputs about surface
% reflectance, illumination, cone sensitivities and monitor phosphors.
%
% surfaces: mxn matrix with each row different surface reflectance function.
% light: 1xn vector of illuminant spectral energy distribution
% cones: 3xn matrix with 3 rows of L,M & S as spectral sensitivies.
% phosphors: 3xn matrix with 3 rows of r,g & b gun spectral
%     energy distributions.
% range: scaling range of color table.  Default is 155.
% header: title of image.  
%
% SEE ALSO:	displayColors.m

% 6/7/96	gmb		Wrote it.

if nargin<5
	range=155;
end

%% Spectral signal from light and surfaces
spectral_signals = surfaces*diag(light);

%% Cone responses to spectral signals.
cone_signals=cones*spectral_signals';

%% Cone responses from monitor phosphors
monitor_to_cones=cones*phosphors';

%% monitor phosphors from cone signals
cones_to_monitor=inv(monitor_to_cones);

%% monitor phosphors from cone signals from surfaces
monitor_signals=cones_to_monitor*cone_signals;

%% correct for gamma
gamma=0.6;
monitor_signals=monitor_signals.^gamma;

%% show the surfaces
if nargin>5
  displayColors(monitor_signals,range,header);
else
  displayColors(monitor_signals,range);
end
