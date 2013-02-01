function [ ] = PlotVr(ConfigFile)
% PlotVr plot azimuthally averaged radial winds
%

[ Config ] = ReadConfig(ConfigFile);

UndefVal = Config.UndefVal;
AzavgDir = Config.AzavgDir;
PlotDir = Config.PlotDir;
ControlCase = Config.ControlCase;

DoDiff = 1;

% limit view into plot
Zmin = 0;
Zmax = 15; %km
Rmin = 40;
Rmax = 70; %km
Tmin = 120;
Tmax = 140;

Phase = 'CO-SS';

% color/contour limits
if (DoDiff == 1)
  Cmin = -5;
  Cmax = 5;
else
  Cmin = -10;
  Cmax = 10;
end

for icase = 1:length(Config.Cases)
  Case = Config.Cases(icase).Cname;
  CaseTitle = regexprep(Case, '.*_', '');
  fprintf('Plotting Vt for case: %s\n', Case)
  fprintf('\n');
  
  % read in the Vr and coordinate data
  % VT data is (x,y,z,t) where x is radius
  Hfile = sprintf('%s/speed_r_%s.h5', AzavgDir, Case);
  fprintf('Reading HDF5 file: %s\n', Hfile);
  VT = squeeze(hdf5read(Hfile, '/speed_r'));
  R = squeeze(hdf5read(Hfile, '/x_coords')) / 1000; % in km
  Z = squeeze(hdf5read(Hfile, '/z_coords')) / 1000; % in km
  T = squeeze(hdf5read(Hfile, '/t_coords')) / 3600; % in hr

  if (DoDiff == 1)
    Hfile = sprintf('%s/speed_r_%s.h5', AzavgDir, ControlCase);
    fprintf('Reading HDF5 file (control): %s\n', Hfile);
    VT_CNTL = squeeze(hdf5read(Hfile, '/speed_r'));
  end

  % trim off radial and vertical ranges
  R1 = find(R >= Rmin, 1, 'first');
  R2 = find(R <= Rmax, 1, 'last');
  Z1 = find(Z >= Zmin, 1, 'first');
  Z2 = find(Z <= Zmax, 1, 'last');
  T1 = find(T >= Tmin, 1, 'first');
  T2 = find(T <= Tmax, 1, 'last');

  % Build data for the plot
  % Convert undefined values to nan so they are excluded from the plot
  % Take average over time (dimension number 3)
  Rvals = R(R1:R2);
  Zvals = Z(Z1:Z2);
  Pdata = squeeze(VT(R1:R2,Z1:Z2,T1:T2));
  Pdata (Pdata  == UndefVal) = nan;
  Pdata = squeeze(nanmean(Pdata,3));

  if (DoDiff == 1)
    PdataCntl = squeeze(VT_CNTL(R1:R2,Z1:Z2,T1:T2));
    PdataCntl (PdataCntl  == UndefVal) = nan;
    PdataCntl = squeeze(nanmean(PdataCntl,3));
    Pdata = Pdata - PdataCntl;
  end

  % plot
  if (DoDiff == 1)
    Pfile  = sprintf('%s/vr_diff_%s.fig', PlotDir, Case);
    Ptitle = sprintf('%s: %s-CLEAN: Vr (m/s)', Phase, CaseTitle);
  else
    Pfile  = sprintf('%s/vr_%s.fig', PlotDir, Case);
    Ptitle = sprintf('%s: %s: Vr (m/s)', Phase, CaseTitle);
  end
  Xlabel = sprintf('Radius (km)');
  Ylabel = sprintf('Height (km)');
   
  Fig = figure;
  
  
  contourf(Rvals, Zvals, Pdata');
  set(gca,'FontSize', 20);
  colormap(redblue);
  shading flat;
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
