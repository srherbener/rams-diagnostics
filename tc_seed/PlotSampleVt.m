function [ ] = PlotSampleVt(ConfigFile)
% PlotSampleVt function to plot Vt
%

[ Config ] = ReadConfig(ConfigFile);

UndefVal = Config.UndefVal;
AzavgDir = Config.AzavgDir;
PlotDir = Config.PlotDir;

Ptime = 90; % plot time - coincides with 

% limit view into plot
Zmin = 0;
Zmax = 15; %km
Rmin = 0;
Rmax = 250; %km

% color/contour limits
Cmin = 0;
Cmax = 80;

% markers
EW_R1 = 0;
EW_R2 = 40; %km

Case = 'TCS_GN_C2000';

fprintf('Plotting Vt for case: %s\n', Case)
fprintf('\n');

% read in the Vt and coordinate data
% VT data is (x,y,z,t) where x is radius
Hfile = sprintf('%s/speed_t_%s.h5', AzavgDir, Case);
fprintf('Reading HDF5 file: %s\n', Hfile);
VT = squeeze(hdf5read(Hfile, '/speed_t'));
R = squeeze(hdf5read(Hfile, '/x_coords')) / 1000; % in km
Z = squeeze(hdf5read(Hfile, '/z_coords')) / 1000; % in km
T = squeeze(hdf5read(Hfile, '/t_coords')) / 3600; % in hr

% trim off radial and vertical ranges
R1 = find(R >= Rmin, 1, 'first');
R2 = find(R <= Rmax, 1, 'last');
Z1 = find(Z >= Zmin, 1, 'first');
Z2 = find(Z <= Zmax, 1, 'last');

TP = find(T >= Ptime, 1, 'first');

% Build data for the plot
% Convert undefined values to nan so they are excluded from the plot
Rvals = R(R1:R2);
Zvals = Z(Z1:Z2);
Pdata = squeeze(VT(R1:R2,Z1:Z2,TP))'; % note transpose - for plot command
Pdata (Pdata  == UndefVal) = nan;

% Create a patch for marking the eyewall region
PX = [  0;
        0;
       40;
       40 ];
PY = [  0;
       15;
       15;
        0 ];
PC = [ 0;
       0;
       0;
       0 ];
PG = 0.8; % Gray scale

% plot
Pfile  = sprintf('%s/ams_vt_SAMPLE_%s.jpg', PlotDir, Case);
Ptitle = sprintf('Vt (m/s)');
Xlabel = sprintf('Radius (km)');
Ylabel = sprintf('Height (km)');
Fsize = 20;
 
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

contourf(Rvals, Zvals, Pdata);
shading flat;
set(gca,'FontSize', Fsize);
colormap(cmap);
caxis([ Cmin Cmax ]);
title(Ptitle);
xlabel(Xlabel);
ylabel(Ylabel);
colorbar;

hold on;

% draw the marker as a "patch"
%patch(PX, PY, PC, 'FaceColor', [ PG PG PG ], 'FaceAlpha', 0.1, 'LineStyle', '--', 'EdgeColor', 'g');
%text(20,13,'EW', 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom', 'FontSize', Fsize);
line([40, 40], [0 15], 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);
text(20,13,'\leftarrow EW \rightarrow', 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom', 'FontSize', 16, 'Color', 'r');

fprintf('Writing plot file: %s\n', Pfile);
saveas(Fig, Pfile);
hold off;
close(Fig);

end
