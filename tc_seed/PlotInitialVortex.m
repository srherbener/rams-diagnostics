function [ ] = PlotInitialVortex(ConfigFile)
% PlotInitialVortex function to plot Vt of initial vortex
%
% This actuall plots Vt at one hour into the simulation since this is the first
% anlysis file available to do this.

[ Config ] = ReadConfig(ConfigFile);

UndefVal = Config.UndefVal;
ControlCase = Config.ControlCase;
AzavgDir = Config.AzavgDir;
PlotDir = Config.PlotDir;

% read in the Vt and coordinate data
Hfile = sprintf('%s/speed_t_%s.h5', AzavgDir, ControlCase)
VT = squeeze(hdf5read(Hfile, '/speed_t'));
X = squeeze(hdf5read(Hfile, '/x_coords'));
Z = squeeze(hdf5read(Hfile, '/z_coords'));

% put X (radius) in km
X = X ./ 1000;

% grab VT at the second time step (first time step is time zero which
% has all wind set to zero).
VT_INIT = squeeze(VT(:,:,2));
VT_INIT = VT_INIT';  % to get radius on x-axis and height on y-axis of
                     % the contour plot

% Convert the undef values to nans so they will be excluded from the plot
VT_INIT(VT_INIT == UndefVal) = nan;

% trim off vertical range since initial vortex Only goes to ~8km
Z1 = find(Z >= 0, 1);
Z2 = find(Z <= 8000, 1, 'last');

Z = Z(Z1:Z2);
VT_INIT = VT_INIT(Z1:Z2,:); % rows are z dimension now (after the transpose)

% plot
Pfile = sprintf('%s/InitVortex.jpg', PlotDir);
Ptitle = sprintf('Initial Vortex: Azimuthally averaged Vt (m/s)');
Xlabel = sprintf('Radius (km)');
Ylabel = sprintf('Height (m)');

Clevs = (0:1:15);

Fig = figure;

contourf(X, Z, VT_INIT, Clevs);
set(gca,'FontSize', 20);
title(Ptitle);
xlabel(Xlabel);
ylabel(Ylabel);
colorbar;

saveas(Fig, Pfile);
close(Fig);

end
