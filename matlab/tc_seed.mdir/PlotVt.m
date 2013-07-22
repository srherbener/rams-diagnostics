function [ ] = PlotVt(ConfigFile)
% PlotVt function to plot Vt of initial vortex
%
% This actuall plots Vt at one hour into the simulation since this is the first
% anlysis file available to do this.

[ Config ] = ReadConfig(ConfigFile);

UndefVal = Config.UndefVal;
ControlCase = Config.ControlCase;
AzavgDir = Config.AzavgDir;
PlotDir = Config.PlotDir;

Stime = 24; % experiment start time
Ptime = 120; % plot time - far enough out to see the size change

% limit view into plot
Zmin = 0;
Zmax = 15; %km
Rmin = 0;
Rmax = 250; %km

% color/contour limits
Cmin = 0;
Cmax = 80;

for icase = 1:length(Config.Cases)
  Case = Config.Cases(icase).Cname;
  CaseTitle = regexprep(Case, '_', '-');
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

  TS = find(T >= Stime, 1, 'first');
  TP = find(T >= Ptime, 1, 'first');

  % Build data for the plot
  % Convert undefined values to nan so they are excluded from the plot
  Rvals = R(R1:R2);
  Zvals = Z(Z1:Z2);
  Pdata = squeeze(VT(R1:R2,Z1:Z2,TP))'; % note transpose - for plot command
  Pdata (Pdata  == UndefVal) = nan;

  % if this is the control case, grab the Vt at the start of the sim too
  if (strcmp(Case,ControlCase))
    PdataStart = squeeze(VT(R1:R2,Z1:Z2,TS))'; % note transpose - for plot command
    PdataStart (PdataStart  == UndefVal) = nan;
  end

  % plot
  Pfile  = sprintf('%s/vt_%s.jpg', PlotDir, Case);
  Ptitle = sprintf('%s: Vt (m/s), T = %d hr', CaseTitle, Ptime);
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
  
  contourf(Rvals, Zvals, Pdata);
  set(gca,'FontSize', 20);
  %colormap(cmap);
  caxis([ Cmin Cmax ]);
  title(Ptitle);
  xlabel(Xlabel);
  ylabel(Ylabel);
  colorbar;
  
  fprintf('Writing plot file: %s\n', Pfile);
  saveas(Fig, Pfile);
  close(Fig);

  if (strcmp(Case,ControlCase))
    Fig = figure;
  
    Pfile  = sprintf('%s/vt_START.jpg', PlotDir);
    Ptitle = sprintf('%s: Vt (m/s), T = %d hr', CaseTitle, Stime);

    contourf(Rvals, Zvals, PdataStart);
    set(gca,'FontSize', 20);
    %colormap(cmap);
    caxis([ Cmin Cmax ]);
    title(Ptitle);
    xlabel(Xlabel);
    ylabel(Ylabel);
    colorbar;
  
    fprintf('Writing plot file: %s\n', Pfile);
    saveas(Fig, Pfile);
    close(Fig);
  end
end

end
