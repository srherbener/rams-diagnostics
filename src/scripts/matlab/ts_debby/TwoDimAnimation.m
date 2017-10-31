function [ ] = TwoDimAnimation(InFile, InVar, OutFile)
% function to create a GIF animation of the contents of a 2D extracted variable

  TWP_DS = ncgeodataset(InFile);
  TWP_VAR = TWP_DS.geovariable(InVar);
  LON_VAR = TWP_DS.geovariable('x_coords');
  LAT_VAR = TWP_DS.geovariable('y_coords');

  TWP = squeeze(TWP_VAR.data(1,:,:));
  LON = LON_VAR.data(:);
  LAT = LAT_VAR.data(:);

  S = TWP_VAR.size;
  Nt = S(1);
  Ny = S(2);
  Nx = S(3);

  if (strcmp(InVar, 'vertint_cond_nd'))
    Clevs = 0.5:0.5:10.5;
    Clims = [ 0.1 5.0 ];
  elseif (strcmp(InVar, 'pi'))
    % Grid 1
    Clevs = 970:2:1010;
    Clims = [ 970 1010 ];

    % Grids 2, 3
    %Clevs = 1000:1:1010;
    %Clims = [ 1000 1010 ];
  elseif (strcmp(InVar, 'theta'))
    % Grid 1 
    Clevs = 285:2:315;
    Clims = [ 285 315 ];

    % Grids 2, 3
    %Clevs = 295:1:305;
    %Clims = [ 290 305 ];
  elseif (strcmp(InVar, 'vapor'))
    Clevs = 0:1:20;
    Clims = [ 0 20 ];
  end

  DelayTime = 0.1;
  Fsize = 22;

  % for Grid3
  LatBounds = [ 7 24 ];
  LonBounds = [ -40 -14 ];

  CoastColor = str2rgb('Black');

  % first establish a graphics context (axis, etc)
  % This plot (outside the loop) will be thrown away, but
  % it is needed to create the "map" variable for the rgb2ind
  % calls inside the loop.
  Fig = figure;
  set(gca,'FontSize', Fsize);

  hold on;
  m_proj('miller', 'lat', LatBounds, 'long', LonBounds);
  m_coast('color', CoastColor, 'linestyle', '-', 'linewidth', 3);
  m_grid('linestyle','none','box','fancy','tickdir','out');
  m_contourf(LON, LAT, TWP, Clevs, 'LineStyle', 'none');

  caxis(Clims);
  colorbar

  set(gca,'nextplot', 'replacechildren');
  f = getframe(gcf);
  [im,map] = rgb2ind(f.cdata,256,'nodither');
  im(1,1,1,Nt) = 0;

  close(Fig);

  % capture frames for the animation
  for it = 1:Nt
    Fig = figure;
    set(gca,'FontSize', Fsize);

    TWP = squeeze(TWP_VAR.data(it,:,:));

    hold on;
    m_proj('miller', 'lat', LatBounds, 'long', LonBounds);
    m_coast('color', CoastColor, 'linestyle', '-', 'linewidth', 3);
    m_grid('linestyle','none','box','fancy','tickdir','out');
    m_contourf(LON,LAT,TWP,Clevs,'LineStyle','none');
    
    colorbar;
    caxis(Clims);

    f = getframe(gcf);
    im(:,:,1,it) = rgb2ind(f.cdata, map, 'nodither');

    close(Fig);
  end  

  imwrite(im, map, OutFile, 'DelayTime', DelayTime, 'LoopCount', 0);
end
