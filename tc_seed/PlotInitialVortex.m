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
Hfile = sprintf('%s/speed_t300_%s.h5', AzavgDir, ControlCase);
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
Z2 = find(Z <= 10000, 1, 'last');

Z = Z(Z1:Z2);
VT_INIT = VT_INIT(Z1:Z2,:); % rows are z dimension now (after the transpose)

% plot
Pfile = sprintf('%s/InitVortex.jpg', PlotDir);
Ptitle = sprintf('Initial Vortex: Azimuthally averaged Vt (m/s)');
Xlabel = sprintf('Radius (km)');
Ylabel = sprintf('Height (m)');

Clevs = (-1:2:21);

Fig = figure;

% create a colormap that goes from dark gray (not black) to white.
%   R == G == B will get a gray shade
%   (R,G,B) = (0,0,0) is black
%   (R,G,B) = (1,1,1) is white
%   need 64 (R,G,B) entries
cband = zeros(64,1);
Cstart = 1;
Cend = 0.40;
Cinc = (Cstart - Cend) / 63;
for i = 1:64
  cband(i,1) = Cstart - ((i-1)*Cinc);
end
cmap = horzcat(cband, cband, cband);
colormap(cmap);


contourf(X, Z, VT_INIT, Clevs);
set(gca,'FontSize', 20);
title(Ptitle);
xlabel(Xlabel);
ylabel(Ylabel);
colorbar;

saveas(Fig, Pfile);
close(Fig);

end
