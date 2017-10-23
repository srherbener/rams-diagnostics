function [ ] = GenTrackLine()

  Pdir = 'Plots';
  if (exist(Pdir, 'dir') ~= 7)
      mkdir(Pdir);
  end
  
  Cases = {
   'TSD_SAL_DUST'
   'TSD_SAL_NODUST'
   'TSD_NONSAL_DUST'
   'TSD_NONSAL_NODUST'
   };
  
   Nc = length(Cases);
  
  % for the TS Debby simulations:
  LatBounds = [ 5 26 ];
  LonBounds = [ -42 -8 ];
  
  % read in the storm location data
  % place the x (Lon) values in a column vector and corresponding
  % y (Lat) values in a separate column vector
  for icase = 1:Nc
    Case = Cases{icase};
    Hfile = sprintf('HDF5/%s/HDF5/storm_center-%s-AS-2006-08-20-120000-g3.h5', Case, Case);
    HdsetLon = '/press_cent_xloc';
    HdsetLat = '/press_cent_yloc';
    fprintf('Reading: %s\n', Hfile);
    fprintf('  Track longitude: %s\n', HdsetLon);
    fprintf('  Track latitude: %s\n', HdsetLat);
  
    % Slons and Slats are (t), column vectors
    Slons = squeeze(h5read(Hfile, HdsetLon));
    Slats = squeeze(h5read(Hfile, HdsetLat));

    if (icase == 1)
      LON = Slons;
      LAT = Slats;
    else
      LON = [ LON' Slons' ]'; % keep as column vectors
      LAT = [ LAT' Slats' ]';
    end
  end
  fprintf('\n');

  % Use polyfit() to do a linear regression
  %   LON is independent variable
  %   LAT is dependent variable
  %
  % After polyfit runs, B(1) is slope, B(2) is y-intercept
  B = polyfit(LON, LAT, 1);

  % evaluate the linear model and generate the R^2 value:
  %   R^2 = 1 - (SSredid/SStotal)
  %    SSresid is sum of squared residuals
  %    SStotal is n-1 * variance (sum of squared differences with mean)
  LAT_FIT = polyval(B, LON);
  SSresid = sum((LAT-LAT_FIT).^2);
  SStotal = (length(LAT)-1) * var(LAT);
  Rsqr = 1 - (SSresid/SStotal);

  fprintf('Results of linear regression:\n');
  fprintf('   Slope: %f\n', B(1));
  fprintf('   Y-Intercept: %f\n', B(2));
  fprintf('\n');
  fprintf('   R-squared: %f\n', Rsqr);
  fprintf('\n');

  % plot
  Fig = figure;
  
  Fsize = 22;
  LegendFsize = 15;

  LineW = 2;
  
  hold on;
  set(gca, 'FontSize', Fsize);

  scatter(LON, LAT, 'o');
  plot(LON, LAT_FIT, 'Color', 'r', 'LineWidth', LineW);

  hold off;
  
%  OutFile = sprintf('%s/TsDebbyTracksAll.jpg', Pdir);
%  fprintf('Writing: %s\n', OutFile);
%  saveas(Fig, OutFile);
%  close(Fig);
end
