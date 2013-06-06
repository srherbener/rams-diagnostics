function [ ] = PlotInitialVortex(ConfigFile)
% PlotInitialVortex function to plot Vt of initial vortex
%
% This actuall plots Vt at one hour into the simulation since this is the first
% anlysis file available to do this.

[ Config ] = ReadConfig(ConfigFile);

UndefVal = Config.UndefVal;
SpinUpCase = Config.SpinUpCase;
AzavgDir = Config.AzavgDir;
PlotDir = Config.PlotDir;

% read in the Vt and coordinate data
%Hfile = sprintf('%s/speed_t300_%s.h5', AzavgDir, SpinUpCase);
Hfile = sprintf('%s/speed_t_%s.h5', AzavgDir, SpinUpCase);
VT = squeeze(hdf5read(Hfile, '/speed_t'));
X = squeeze(hdf5read(Hfile, '/x_coords'))/1000; % in km
Z = squeeze(hdf5read(Hfile, '/z_coords'))/1000; % in km

% grab VT at the second time step (first time step is time zero which
% has all wind set to zero).
VT_INIT = squeeze(VT(:,:,2));
VT_INIT = VT_INIT';  % to get radius on x-axis and height on y-axis of
                     % the contour plot

% Convert the undef values to nans so they will be excluded from the plot
VT_INIT(VT_INIT == UndefVal) = nan;

% trim off vertical range since initial vortex Only goes to ~8km
Z1 = find(Z >= 0, 1);
Z2 = find(Z <= 10, 1, 'last');

Z = Z(Z1:Z2);
VT_INIT = VT_INIT(Z1:Z2,:); % rows are z dimension now (after the transpose)

% plot
Pfile = sprintf('%s/InitVortex.jpg', PlotDir);
Xlabel = sprintf('Radius (km)');
Ylabel = sprintf('Height (km)');

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


contourf(X, Z, VT_INIT);
shading flat;
set(gca,'FontSize', 45);
caxis([ 0 15 ]);
set(gca, 'XTick', [ 100 200 ]);
set(gca, 'LineWidth', 2);
set(gca, 'TickLength', [ 0.025 0.025 ]);
% The title is in a box that adjusts to the amount of characters in
% the title. Ie, it doesn't do any good to do Left/Center/Right
% alignment. But, the entire box can be moved to the left side of the
% plot.
T = title('(a)');
set(T, 'Units', 'Normalized');
set(T, 'HorizontalAlignment', 'Left');
Tpos = get(T, 'Position');
Tpos(1) = 0; % line up with left edge of plot area
set(T, 'Position', Tpos);
xlabel(Xlabel);
ylabel(Ylabel);
colorbar;

% Fix up the positioning
Ppos = get(gca, 'Position'); % position of plot area
Ppos(1) = Ppos(1) * 1.00;
Ppos(2) = Ppos(2) * 0.90;
Ppos(3) = Ppos(3) * 0.90;
Ppos(4) = Ppos(4) * 0.90;
set(gca, 'Position', Ppos);

saveas(Fig, Pfile);
close(Fig);

end
