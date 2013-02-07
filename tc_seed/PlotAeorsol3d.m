function [ ] = PlotAeorsol3d(ConfigFile, ReadRevu)
% PlotAeorsol3d function to make 3d plot of TC ingesting aerosols
%
%  ReadRevu: 1 --> read in and select data from REVU files, then
%                  store in intermediate file.
%            0 --> assume intermediate file is built and start by
%                  reading data from the intermediate file.
%                   

  [ Config ] = ReadConfig(ConfigFile);
  
  Pdir = Config.PlotDir;
  Ddir = Config.DiagDir;
  
  % Make sure output directories exists
  if (exist(Pdir, 'dir') ~= 7)
      mkdir(Pdir);
  end
  
  if (exist(Ddir, 'dir') ~= 7)
      mkdir(Ddir);
  end
  
  % big font
  Fsize = 20;

  
  % Ranges for selection of data
  %   Domain center is 15N, 40W (longitude is ref to East, so this is -40)
  Xmin = -42;
  Xmax = -38;
  Ymin = 13;
  Ymax = 16;
  Zmin = 0;
  Zmax = 6500;
  Tplot = 100 * 3600;

  % Iso surface values
  CloudIsurf   = 0.25;
  CloudAlpha   = 0.3;
  AerosolIsurf = 200;
  AerosolAlpha = 0.8;
  CcnHydroIsurf = 0.5;
  CcnHydroAlpha = 0.7;
  
  IngestView = [ -60 90 ];
  ScavengeView = [ - 60 90 ];

  % Use 2000/cc case for now - REVU output files
  % Want to plot cloud mass and aerosol number
  CloudFile   = 'cloud-TCS_SD_C2000-AS-1998-08-22-120000-g3.h5';
  CloudDset   = 'cloud';
  AerosolFile = 'ccn_conc-TCS_SD_C2000-AS-1998-08-22-120000-g3.h5';
  AerosolDset = 'ccn_concen';
  CcnHydroFile = 'ccn_rain_mass-TCS_SD_C2000-AS-1998-08-22-120000-g3.h5';
  CcnHydroDset = 'ccn_rain_mass';
  
  % Intermediate file
  DataFile = sprintf('%s/Aerosol3d.h5', Ddir);
  CldDset  = 'Cloud';
  AeroDset = 'Aerosol';
  LatDset  = 'Lat';
  LonDset  = 'Lon';
  HgtDset  = 'Height';

  fprintf('***********************************************************************\n');
  fprintf('Plotting Aerosol Ingestion:\n');
  fprintf('\n');
  
  % Only build the intermediate file if ReadRevu is a one. This way
  % we can skip the long read time from REVU files when fiddling with
  % the plots.
  if (ReadRevu == 1)
    % Read in the data one var at a time since the datasets can get big.
    fprintf('Building data file from REVU files:\n');
    fprintf('  Reading file: %s, Dataset: %s\n', CloudFile, CloudDset);
    CLOUD  = ReadSelectXyzt(CloudFile, CloudDset, Xmin, Xmax, Ymin, Ymax, Zmin, Zmax, Tplot, Tplot);
    LON    = ReadSelectXyzt(CloudFile, 'x_coords', Xmin, Xmax, Ymin, Ymax, Zmin, Zmax, Tplot, Tplot);
    LAT    = ReadSelectXyzt(CloudFile, 'y_coords', Xmin, Xmax, Ymin, Ymax, Zmin, Zmax, Tplot, Tplot);
    HGT    = ReadSelectXyzt(CloudFile, 'z_coords', Xmin, Xmax, Ymin, Ymax, Zmin, Zmax, Tplot, Tplot);

    fprintf('  Reading file: %s, Dataset: %s\n', AerosolFile, AerosolDset);
    AERO = ReadSelectXyzt(AerosolFile, AerosolDset, Xmin, Xmax, Ymin, Ymax, Zmin, Zmax, Tplot, Tplot);

    fprintf('  Reading file: %s, Dataset: %s\n', CcnHydroFile, CcnHydroDset);
    CCN_HYDRO = ReadSelectXyzt(CcnHydroFile, CcnHydroDset, Xmin, Xmax, Ymin, Ymax, Zmin, Zmax, Tplot, Tplot);

    % Write out to the intermediate data file
    fprintf('  Writing file: %s\n', DataFile);
    hdf5write(DataFile, CldDset, CLOUD);
    hdf5write(DataFile, AeroDset, AERO, 'WriteMode', 'append');
    hdf5write(DataFile, CcnHydroDset, CCN_HYDRO, 'WriteMode', 'append');
    hdf5write(DataFile, LonDset, LON, 'WriteMode', 'append');
    hdf5write(DataFile, LatDset', LAT, 'WriteMode', 'append');
    hdf5write(DataFile, HgtDset, HGT, 'WriteMode', 'append');
    fprintf('\n');
  end

  % read in from intermediate file and do the plot
  fprintf('Reading data file: %s\n', DataFile);
  fprintf('\n');
  CLOUD     = hdf5read(DataFile, CldDset);
  AERO      = hdf5read(DataFile, AeroDset);
  CCN_HYDRO = hdf5read(DataFile, CcnHydroDset);
  LON       = hdf5read(DataFile, LonDset);
  LAT       = hdf5read(DataFile, LatDset);
  HGT       = hdf5read(DataFile, HgtDset)/1000; % change to km

  % for plotting
  CLOUD = permute(CLOUD, [ 2 1 3 ]);
  AERO = permute(AERO, [ 2 1 3 ]);
  CCN_HYDRO = permute(CCN_HYDRO, [ 2 1 3 ]);
  
  % Plot ingesting
  OutFile = sprintf('%s/TcAerosolIngest.fig', Pdir);
  Ptitle = sprintf('TC Ingesting Aerosols');
  Xlabel = 'Lon (deg)';
  Ylabel = 'Lat (deg)';
  Zlabel = 'Height (km)';
  
  % Create a grid
  [ LONG, LATG, HGTG ] = meshgrid(LON, LAT, HGT);
  
  % make the cloud data transparent so you can see the 
  % aerosol underneath
  Fig = figure;
  PC = patch(isosurface(LONG, LATG, HGTG, CLOUD, CloudIsurf));
  set(gca, 'FontSize', Fsize);
  isonormals(LONG, LATG, HGTG, CLOUD, PC);
  set(PC,'FaceColor', 'Cyan', 'EdgeColor', 'None', 'FaceAlpha', CloudAlpha);
  hold on;

  PA = patch(isosurface(LONG, LATG, HGTG, AERO, AerosolIsurf));
  isonormals(LONG, LATG, HGTG, AERO, PA);
  set(PA,'FaceColor', 'Red', 'EdgeColor', 'None', 'FaceAlpha', AerosolAlpha);

  grid on;
  view(3);
  %view(gca, [-61.5 22]); % changes rotation and height, [ azimuth, elevation ]
  view(gca, IngestView); % changes rotation and height
  camlight;
  lighting gouraud;
  
  title(Ptitle);
  xlabel(Xlabel);
  ylabel(Ylabel);
  zlabel(Zlabel);

  fprintf('Writing file: %s\n', OutFile);
  fprintf('\n');
  saveas(Fig, OutFile);
  close(Fig);

  
  % Plot scavenging
  Ptitle = sprintf('TC Scavenging out Aerosols');
  OutFile = sprintf('%s/TcAerosolScavenge.fig', Pdir);
  Xlabel = 'Lon (deg)';
  Ylabel = 'Lat (deg)';
  Zlabel = 'Height (km)';
  
  % Create a grid
  [ LONG, LATG, HGTG ] = meshgrid(LON, LAT, HGT);
  
  % make the cloud data transparent so you can see the 
  % aerosol underneath
  Fig = figure;
  PC = patch(isosurface(LONG, LATG, HGTG, CLOUD, CloudIsurf));
  set(gca, 'FontSize', Fsize);
  isonormals(LONG, LATG, HGTG, CLOUD, PC);
  set(PC,'FaceColor', 'Cyan', 'EdgeColor', 'None', 'FaceAlpha', CloudAlpha);
  hold on;

  PA = patch(isosurface(LONG, LATG, HGTG, CCN_HYDRO, CcnHydroIsurf));
  isonormals(LONG, LATG, HGTG, CCN_HYDRO, PA);
  set(PA,'FaceColor', 'Red', 'EdgeColor', 'None', 'FaceAlpha', CcnHydroAlpha);

  grid on;
  view(3);
  %view(gca, [-61.5 22]); % changes rotation and height, [ azimuth, elevation ]
  view(gca, ScavengeView); % changes rotation and height
  camlight;
  lighting gouraud;
  
  title(Ptitle);
  xlabel(Xlabel);
  ylabel(Ylabel);
  zlabel(Zlabel);

  fprintf('Writing file: %s\n', OutFile);
  fprintf('\n');
  saveas(Fig, OutFile);
  close(Fig);
end
