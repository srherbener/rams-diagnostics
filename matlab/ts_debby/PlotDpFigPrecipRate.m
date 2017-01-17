function [ ] = PlotDpFigPrecipRate()

  Pdir = 'Plots';
  if (exist(Pdir, 'dir') ~= 7)
      mkdir(Pdir);
  end
  
  Fsize = 25;

  % Storm center data.
  InFile = 'HDF5/TSD_SAL_DUST/HDF5/storm_center-TSD_SAL_DUST-AS-2006-08-20-120000-g3.h5';

  % Read in storm center longitudes and latitudes
  InVname = '/press_cent_xloc';
  fprintf('Reading: %s (%s)\n', InFile, InVname);
  TSD_LON_LOC = double(squeeze(h5read(InFile, InVname)));

  InVname = '/press_cent_yloc';
  fprintf('Reading: %s (%s)\n', InFile, InVname);
  TSD_LAT_LOC = double(squeeze(h5read(InFile, InVname)));

  % Read in the storm location for the following times
  %     t =  31 -> 21Z, Aug22
  %     t =  61 -> 12Z, Aug23
  %     t =  91 -> 03Z, Aug24
  %     t = 121 -> 18Z, Aug24
  TSD_LOC1 = [ TSD_LAT_LOC( 31) TSD_LON_LOC( 31) ];
  TSD_LOC2 = [ TSD_LAT_LOC( 61) TSD_LON_LOC( 61) ];
  TSD_LOC3 = [ TSD_LAT_LOC( 91) TSD_LON_LOC( 91) ];
  TSD_LOC4 = [ TSD_LAT_LOC(121) TSD_LON_LOC(121) ];

  Ptitle1 = '21Z, 22Aug';
  Ptitle2 = '12Z, 23Aug';
  Ptitle3 = '03Z, 24Aug';
  Ptitle4 = '24Z, 24Aug';

  % Get grid 3 location
  GlocFname = 'GridLocs.txt';
  Grid3Loc = ReadGridLoc(GlocFname, 'Grid3:');

  % Map locations:
  %  grid3 extended a small amount
  MapLoc = Grid3Loc + [ -2 +3 -2 +6 ];

  % Dust deposited on the surface
  InFile = 'HDF5/TSD_SAL_DUST/HDF5/pcprate-TSD_SAL_DUST-AS-2006-08-20-120000-g3.h5';
  InVname = '/pcprate';

  % Read in and create snapshots of the surface deposited dust
  fprintf('Reading: %s (%s)\n', InFile, InVname);
  LON    = double(squeeze(h5read(InFile, '/x_coords')));
  LAT    = double(squeeze(h5read(InFile, '/y_coords')));
  Nlon = length(LON);
  Nlat = length(LAT);
  % pcprate is (x,y,t)
  %   times
  %     t =  31 -> 21Z, Aug22
  %     t =  61 -> 12Z, Aug23
  %     t =  91 -> 03Z, Aug24
  %     t = 121 -> 18Z, Aug24
  %   units are mm/h
  %
  % Note transpose of data (in order to get LAT, LON, PRATE to line up for
  % the contourfm command.
  PRATE1 = double(squeeze(h5read(InFile, InVname,[ 1 1  31 ], [ Nlon Nlat 1 ])))';
  PRATE2 = double(squeeze(h5read(InFile, InVname,[ 1 1  61 ], [ Nlon Nlat 1 ])))';
  PRATE3 = double(squeeze(h5read(InFile, InVname,[ 1 1  91 ], [ Nlon Nlat 1 ])))';
  PRATE4 = double(squeeze(h5read(InFile, InVname,[ 1 1 121 ], [ Nlon Nlat 1 ])))';

  % Plot

  % Contour levels
  Cmin = 0;
  Cinc = 0.1;
  Cmax = 1;

  Cmap = 'parula';

  % Create each panel in a separate figure since can only have one axesm per figure
  [ F1 F1axes ] = PlotDpFigMapContour(MapLoc, LAT, LON, PRATE1, TSD_LOC1, 'a', Ptitle1, Cmin, Cinc, Cmax, Cmap, Fsize);
  [ F2 F2axes ] = PlotDpFigMapContour(MapLoc, LAT, LON, PRATE2, TSD_LOC2, 'b', Ptitle2, Cmin, Cinc, Cmax, Cmap, Fsize);
  [ F3 F3axes ] = PlotDpFigMapContour(MapLoc, LAT, LON, PRATE3, TSD_LOC3, 'c', Ptitle3, Cmin, Cinc, Cmax, Cmap, Fsize);
  [ F4 F4axes ] = PlotDpFigMapContour(MapLoc, LAT, LON, PRATE4, TSD_LOC4, 'd', Ptitle4, Cmin, Cinc, Cmax, Cmap, Fsize);

  % Save each panel in a separate file, for now
  OutFile = sprintf('%s/DpFigPrecipRate1.jpg', Pdir);
  fprintf('Writing: %s\n', OutFile);
  saveas(F1, OutFile);
  close(F1);

  OutFile = sprintf('%s/DpFigPrecipRate2.jpg', Pdir);
  fprintf('Writing: %s\n', OutFile);
  saveas(F2, OutFile);
  close(F2);

  OutFile = sprintf('%s/DpFigPrecipRate3.jpg', Pdir);
  fprintf('Writing: %s\n', OutFile);
  saveas(F3, OutFile);
  close(F3);

  OutFile = sprintf('%s/DpFigPrecipRate4.jpg', Pdir);
  fprintf('Writing: %s\n', OutFile);
  saveas(F4, OutFile);
  close(F4);

end
