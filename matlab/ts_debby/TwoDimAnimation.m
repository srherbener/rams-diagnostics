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
    Clevs = 1000:1:1010;
    Clims = [ 1000 1010 ];
  elseif (strcmp(InVar, 'theta'))
    Clevs = 295:1:305;
    Clims = [ 290 305 ];
  elseif (strcmp(InVar, 'vapor'))
    Clevs = 0:1:20;
    Clims = [ 0 20 ];
  end

  DelayTime = 0.1;
  Fsize = 25;

  % first establish a graphics context (axis, etc)
  contourf(LON, LAT, TWP, Clevs, 'LineStyle', 'none');
  shading flat;
%  axis tight;
%  set(gca,'nextplot','replacechildren','visible','off');
  set(gca,'FontSize', Fsize);
  title('TS Debby (2006)');
  xlabel('Lon'); 
  ylabel('Lat'); 
  caxis(Clims);
  colorbar

  set(gca,'nextplot', 'replacechildren');
  f = getframe(gcf);
  [im,map] = rgb2ind(f.cdata,256,'nodither');
  im(1,1,1,Nt) = 0;

  % capture frames for the animation
  for it = 1:Nt
    TWP = squeeze(TWP_VAR.data(it,:,:));
    contourf(LON, LAT, TWP, Clevs, 'LineStyle', 'none');
    shading flat;
    colorbar;
    caxis(Clims);
    f = getframe(gcf);
    im(:,:,1,it) = rgb2ind(f.cdata, map, 'nodither');
  end  

  imwrite(im, map, OutFile, 'DelayTime', DelayTime, 'LoopCount', 0);
end
