function [ ] = PlotDpFigDustSfcDep()

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
  InFile = 'HDF5/TSD_SAL_DUST/HDF5/dust_sfc-TSD_SAL_DUST-AS-2006-08-20-120000-g3.h5';
  InVname = '/accpdust';

  % Read in and create snapshots of the surface deposited dust
  fprintf('Reading: %s (%s)\n', InFile, InVname);
  LON    = double(squeeze(h5read(InFile, '/x_coords')));
  LAT    = double(squeeze(h5read(InFile, '/y_coords')));
  Nlon = length(LON);
  Nlat = length(LAT);
  % accpdust is (x,y,t)
  %   times
  %     t =  31 -> 21Z, Aug22
  %     t =  61 -> 12Z, Aug23
  %     t =  91 -> 03Z, Aug24
  %     t = 121 -> 18Z, Aug24
  %   units are mg/m^2
  %   convert to kg/km^2
  %      accpdust in mg/m^2 * 10e6 m^2/km^2 * kg/10e6 mg -> accpdust in kg/km^2
  %
  % Note transpose of data (in order to get LAT, LON, MD_SFC to line up for
  % the contourfm command.
  MD_SFC1 = double(squeeze(h5read(InFile, InVname,[ 1 1  31 ], [ Nlon Nlat 1 ])))';
  MD_SFC2 = double(squeeze(h5read(InFile, InVname,[ 1 1  61 ], [ Nlon Nlat 1 ])))';
  MD_SFC3 = double(squeeze(h5read(InFile, InVname,[ 1 1  91 ], [ Nlon Nlat 1 ])))';
  MD_SFC4 = double(squeeze(h5read(InFile, InVname,[ 1 1 121 ], [ Nlon Nlat 1 ])))';

  % Plot

  % Create each panel in a separate figure since can only have one axesm per figure
  [ F1 F1axes ] = PlotDpFigMapContour(MapLoc, LAT, LON, MD_SFC1, TSD_LOC1, 'a', Ptitle1, Fsize);
  [ F2 F2axes ] = PlotDpFigMapContour(MapLoc, LAT, LON, MD_SFC2, TSD_LOC2, 'b', Ptitle2, Fsize);
  [ F3 F3axes ] = PlotDpFigMapContour(MapLoc, LAT, LON, MD_SFC3, TSD_LOC3, 'c', Ptitle3, Fsize);
  [ F4 F4axes ] = PlotDpFigMapContour(MapLoc, LAT, LON, MD_SFC4, TSD_LOC4, 'd', Ptitle4, Fsize);

  % Save each panel in a separate file, for now
  OutFile = sprintf('%s/DpFigDustSfcDep1.jpg', Pdir);
  fprintf('Writing: %s\n', OutFile);
  saveas(F1, OutFile);
  close(F1);

  OutFile = sprintf('%s/DpFigDustSfcDep2.jpg', Pdir);
  fprintf('Writing: %s\n', OutFile);
  saveas(F2, OutFile);
  close(F2);

  OutFile = sprintf('%s/DpFigDustSfcDep3.jpg', Pdir);
  fprintf('Writing: %s\n', OutFile);
  saveas(F3, OutFile);
  close(F3);

  OutFile = sprintf('%s/DpFigDustSfcDep4.jpg', Pdir);
  fprintf('Writing: %s\n', OutFile);
  saveas(F4, OutFile);
  close(F4);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ Fig Paxes ] = PlotDpFigMapContour(MapLoc, LAT, LON, Z, StormLoc, Pmarker, Ptitle, Fsize)

  Fig = figure;

  % contour levels
  Cmin = 0;
  Cinc = 1;
  Cmax = 10;
  Clevs = [ Cmin:Cinc:Cmax ];

  LatBounds = [ MapLoc(1) MapLoc(2) ];
  LonBounds = [ MapLoc(3) MapLoc(4) ];

  StormLat = [ StormLoc(1) ];
  StormLon = [ StormLoc(2) ];

  LineW = 3;
  Aratio = 1.4;

  % Location on page
  PlotLoc = [ 0.2 0.2 0.6 0.7 ];
  CbarLoc = [ 0.2 0.1 0.6 0.1 ];

  % Use cylindrical projection - Miller
  % Change aspect ratio of plotbox and data to make map longer in the x axes direction
  Paxes = worldmap(LatBounds, LonBounds);
  setm(Paxes, 'MapProjection', 'miller');
  set(Paxes, 'Position', PlotLoc);
  set(Paxes, 'FontSize', Fsize);
  setm(Paxes, 'FontSize', Fsize);
  pbaspect( [ Aratio 1 1 ]);
  daspect( [ 1 Aratio 1 ]);

  % Plot on linear scale
  contourfm(LAT, LON, Z, Clevs, 'LineStyle', 'none'); 
  caxis([ Cmin Cmax ]);

  % Mark TS Debby location
  plotm(StormLat, StormLon, 'Color', 'r', 'LineStyle', 'none', 'Marker', 'x', 'LineWidth', LineW, 'MarkerSize', Fsize);

  % Colorbar, this sets the colormap too.
  Cbar = contourcmap('hot', 'ColorBar', 'on', 'Location', 'horizontal');
  set(Cbar, 'Position', CbarLoc);
  set(Cbar, 'FontSize', Fsize);
  xlabel(Cbar, 'M_D_S_F_C (kg km^-^2)');

%   % Reduce the amount of tick marks so that the labels won't overlap
%   Cticks = [ Cmin:Cinc:Cmax ];
%   CtickLabels = { '0.0' '0.6' '1.2' '1.8' '2.4' '3.0' };
%   set(Cbar, 'XTick', Cticks);
%   set(Cbar, 'XTickLabel', CtickLabels);

  % Reverse the 'hot' colormap for showing dust
  % Make sure this code comes after the contourcmap call since contourcmap
  % also sets the colormap.
  Cmap = colormap('hot');
  Cmap = Cmap(end:-1:1,:); % reverse the entries
  colormap(Paxes, Cmap);
  colormap(Cbar, Cmap);

  % title - left justify if a panel marker was specified
  if (isempty(Pmarker))
    title(Ptitle);
  else
    T = title(sprintf('(%s) %s', Pmarker, Ptitle));
    LeftJustTitle(Paxes, T);
  end
end
