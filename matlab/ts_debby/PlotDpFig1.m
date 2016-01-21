function [ ] = PlotDpFig1()
% PlotDpFig1 

  Pdir = 'Plots';
  if (exist(Pdir, 'dir') ~= 7)
      mkdir(Pdir);
  end
  
  Fsize = 18;
  
  % Images
  SalPic1 = 'IMAGES/SAL_Aug22_06Z.png'; 
  SalPic2 = 'IMAGES/SAL_Aug23_12Z.png'; 
  SalPic3 = 'IMAGES/SAL_Aug24_18Z.png'; 

  % Vertically integrated dust mass
  VintDustFile = 'HDF5/TSD_SAL_DUST/HDF5/vint_dust-TSD_SAL_DUST-AS-2006-08-20-120000-g3.h5';
  VintDustVname = '/vertint_dust';

  T1 = 1;    % sample time step numbers
  T2 = 61; 
  T3 = 121; 

  % Read in and create 3 snapshots of the vertically integrated dust field
  fprintf('Reading: %s (%s)\n', VintDustFile, VintDustVname);
  DUST = squeeze(h5read(VintDustFile, VintDustVname));
  X    = squeeze(h5read(VintDustFile, '/x_coords'));
  Y    = squeeze(h5read(VintDustFile, '/y_coords'));
  Z    = squeeze(h5read(VintDustFile, '/z_coords'));
  T    = squeeze(h5read(VintDustFile, '/t_coords'));

  Nx = length(X);
  Ny = length(Y);
  Nz = length(Z);
  Nt = length(T);

  D1 = log10(squeeze(DUST(:,:,T1)));
  D2 = log10(squeeze(DUST(:,:,T2)));
  D3 = log10(squeeze(DUST(:,:,T3)));

  
  % plot
  OutFile = sprintf('%s/DpFig1_StormSeq.jpg', Pdir);
  
  Fig = figure;

  % Place SAL images into panels on the left side of figure
  Axes = subplot(3,2,1);
  PlaceSalImage(Axes, SalPic1, 'a', 'Aug22:06Z', Fsize);

  Axes = subplot(3,2,3);
  PlaceSalImage(Axes, SalPic2, 'c', 'Aug23:12Z', Fsize);

  Axes = subplot(3,2,5);
  PlaceSalImage(Axes, SalPic3, 'e', 'Aug24:18Z', Fsize);

  % Plot simulations (contourf) in panels on the right side of figure
  Axes = subplot(3,2,2);
  PlotVintDust(Axes, X, Y, D1, 'b', 'Aug22:06Z', Fsize, 0, 1);

  Axes = subplot(3,2,4);
  PlotVintDust(Axes, X, Y, D2, 'd', 'Aug23:12Z', Fsize, 0, 1);

  Axes = subplot(3,2,6);
  PlotVintDust(Axes, X, Y, D3, 'f', 'Aug24:18Z', Fsize, 1, 1);
  
  fprintf('Writing: %s\n', OutFile);
  saveas(Fig, OutFile);
  close(Fig);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = PlaceSalImage(Paxes, ImFile, Pmarker, Ptitle, Fsize)

  axes(Paxes);
  % Axes position is the location of the axes within the entire figure.
  % 0 -> left side, or bottom
  % 1 -> right side, or top
  %  position -> [ left_x bottom_y width_x height_y ]
  Apos = get(Paxes, 'Position');
  ImgScale = 1.2; % make image 12% larger
  Apos(1) = 0.05;                % shift to the left
  Apos(2) = Apos(2) - 0.02;      % shift downward
  Apos(3) = Apos(3) * ImgScale;
  Apos(4) = Apos(4) * ImgScale;

  set(Paxes, 'Position', Apos);
  imshow(ImFile);

  set(Paxes, 'FontSize', Fsize);

  if (isempty(Pmarker))
    title(Ptitle);
  else
    T = title(sprintf('(%s) %s', Pmarker, Ptitle));
    LeftJustTitle(T);
  end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = PlotVintDust(Paxes, X, Y, Z, Pmarker, Ptitle, Fsize, ShowX, ShowY)

  axes(Paxes);

  contourf(X, Y, Z', 'LineStyle', 'none');
  Cbar = colorbar;
  caxis([ -2 6 ]);

  set(Cbar, 'Ticks', [ -1 2 5 ]);
  set(Cbar, 'TickLabels', { '10^-^2' '10^2' '10^5' });

  set(Paxes, 'FontSize', Fsize);
  if (ShowX > 0)
    xlabel('Lon (\circ)');
    set(Paxes, 'XTick', [ -35 -25 -15 ]);
  else
    set(Paxes, 'XTickLabel', {});
  end
  if (ShowY > 0)
    ylabel('Lat (\circ)');
    set(Paxes, 'YTick', [ 10 15 20 ]);
  else
    set(Paxes, 'YTickLabel', {});
  end

  if (isempty(Pmarker))
    title(Ptitle);
  else
    T = title(sprintf('(%s) %s', Pmarker, Ptitle));
    LeftJustTitle(T);
  end

end
