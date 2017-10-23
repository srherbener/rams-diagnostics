function [ ] = GenInitPlots(IcFile, ScFile)
% GenInitPlots generate plots containing initial conditions for Thorncroft et al., 1993 ETC lifecycles

  % Norhtern Hemisphere winter conditions, pick an arbitrary date in winter
  fprintf('Generating Plots with ETC initial conditions:\n');
  fprintf('  Plot file: %s\n', IcFile);
  fprintf('  Plot file: %s\n', ScFile);

  % Input Files
  Ufile = 'HDF5/u-a-AS-2014-01-01-120000-g1.h5';
  Vfile = 'HDF5/v-a-AS-2014-01-01-120000-g1.h5';
  VaporFile = 'HDF5/vapor-a-AS-2014-01-01-120000-g1.h5';
  ThetaFile = 'HDF5/theta-a-AS-2014-01-01-120000-g1.h5';
  SpFile = 'HDF5/sea_press-a-AS-2014-01-01-120000-g1.h5';

  Xmid = 80;
  

  % Read in data for initial conditions plot
  U_DS = ncgeodataset(Ufile);
  V_DS = ncgeodataset(Vfile);
  VAP_DS = ncgeodataset(VaporFile);
  THETA_DS = ncgeodataset(ThetaFile);
  SP_DS = ncgeodataset(SpFile);

  U_VAR = U_DS.geovariable('u');
  V_VAR = V_DS.geovariable('v');
  VAP_VAR = VAP_DS.geovariable('vapor');
  THETA_VAR = THETA_DS.geovariable('theta');
  SP_VAR = SP_DS.geovariable('sea_press');

  % coords from the u-file
  LAT_VAR = U_DS.geovariable('y_coords');
  LON_VAR = U_DS.geovariable('x_coords');
  Z_VAR = U_DS.geovariable('z_coords');

  %%%%%%%%%%%%%%%%% INITIAL CONDITIONS %%%%%%%%%%%%%%%%%%%%
  % data is organized as (t, z, y, x)
  U = squeeze(U_VAR.data(1, :, :, Xmid));
  VAP = squeeze(VAP_VAR.data(1, :, :, Xmid));
  THETA = squeeze(THETA_VAR.data(1, :, :, Xmid));

  LAT = LAT_VAR.data(:);      % Degrees
  LON = LON_VAR.data(:);      % Degrees
  Z = Z_VAR.data(:) ./ 1000;  % km

  % Place u on filled contours
  % Theta on contour lines
  Fig = figure;

  Fsize = 20;
  CbarFsize = 20;
  Cfontsize = 20;
  Clabspace = 500;
  Clinewidth = 1.2;

  Uclevs = -20:5:40;
  ThetaClevs = 250:50:700;
  VapClevs = 0:3:18;

  % U
  ax(1) = axes;
  [C1, h1 ] = contourf(LAT, Z, U, Uclevs);
  set(ax(1), 'FontSize', Fsize);
  shading flat;
  cbar = colorbar;
  set(cbar, 'FontSize', CbarFsize);
  title('(a) Initial Conditions');
  xlabel('Latitude');
  ylabel('Height (km)');
  % the different axes labeling, fontsize, colorbar change the axes position
  % therefore copy the positions of axes 1 to the other axes
  % record the position last in the sequence of plotting commands for U
  Cpos = get(ax(1), 'Position');

  % Theta
  ax(2) = axes; 
  [ C2, h2 ] = contour(LAT, Z, THETA, ThetaClevs, 'Color', 'k', 'LineWidth', Clinewidth);
  set(ax(2), 'Visible', 'off');
  set(ax(2), 'Position', Cpos);
  clabel(C2, h2, ThetaClevs(2:2:end), 'Color', 'k', 'FontSize', Cfontsize, 'LabelSpacing', Clabspace);

  % Vapor
  ax(3) = axes; 
  [ C3, h3 ] = contour(LAT, Z, VAP, VapClevs, 'Color', 'r', 'LineWidth', Clinewidth);
  set(ax(3), 'Visible', 'off');
  set(ax(3), 'Position', Cpos);
  clabel(C3, h3, VapClevs(2:2:end), 'Color', 'r', 'FontSize', Cfontsize, 'LabelSpacing', Clabspace);

  saveas(Fig, IcFile);
  close(Fig);

  %%%%%%%%%%%%%%%%% SYNOPTIC CONDITIONS %%%%%%%%%%%%%%%%%%%%
  % U and V are (t,z,y,x)
  % sea_press is (t,y,x)

  Tend = 121; % 5 days
  Zupper = 41; % z = 12km (~ 200mb)

  WindClevs = 0:5:50;
  SpClevs = 990:5:1020;

  clear ax;
  clear Cpos;
  

  U = squeeze(U_VAR.data(Tend, Zupper, :, :));
  V = squeeze(V_VAR.data(Tend, Zupper, :, :));
  SP = squeeze(SP_VAR.data(Tend,:,:));

  % Form magnitude of wind
  WIND_SPEED = sqrt(U.^2 + V.^2);

  Fig = figure;

  % Upper level wind speed
  ax(1) = axes;
  [ C1, h1 ] = contourf(LON, LAT, WIND_SPEED, WindClevs);
  set(ax(1), 'FontSize', Fsize);
  shading flat;
  cbar = colorbar;
  set(cbar, 'FontSize', CbarFsize);
  title('(b) Synoptic Conditions, Day 5'); 
  xlabel('Longitude');
  ylabel('Latitude');
  % the different axes labeling, fontsize, colorbar change the axes position
  % therefore copy the positions of axes 1 to the other axes
  % record the position last in the sequence of plotting commands for U
  Cpos = get(ax(1), 'Position');

  % Theta
  ax(2) = axes; 
  [ C2, h2 ] = contour(LON, LAT, SP, SpClevs, 'Color', 'c', 'LineWidth', Clinewidth);
  set(ax(2), 'Visible', 'off');
  set(ax(2), 'Position', Cpos);
  clabel(C2, h2, SpClevs(2:2:end), 'Color', 'c', 'FontSize', Cfontsize, 'LabelSpacing', Clabspace);


  saveas(Fig, ScFile);
  close(Fig);
end
