function [] = AeorsolAnimation ( InFile, Var, TsInc, Delay, Zlevel, Pmin, Pmax, Ptitle, OutFile )
% AeorsolAnimation Create animated GIF for aerosol source

  Paspect = [ 1 1 1 ];
  Tstart = 1;
  Tend = 110;

  fprintf('Creating GIF movie:\n');
  fprintf('  Input file: %s\n', InFile);
  fprintf('  Variable: %s\n', Var);
  fprintf('  Time step increment: %d\n', TsInc);
  fprintf('  Movie frame delay: %.2f\n', Delay);
  fprintf('  Plot minimum: %.2f\n', Pmin);
  fprintf('  Plot maximum: %.2f\n', Pmax);
  fprintf('  Output File: %s\n', OutFile);
  fprintf('\n');

  % Read in the variable
  % Data is (x,y,z,t)
  %   x - longitude
  %   y - latitude
  %   z - levels
  %   t - time
  fprintf('Reading HDF5 file: %s, Dataset: %s\n', InFile, Var);
  fprintf('\n');
  Hvar = hdf5read(InFile, Var);
  Lon = hdf5read(InFile, 'x_coords');
  Lat = hdf5read(InFile, 'y_coords');
  Z = hdf5read(InFile, 'z_coords') / 1000; % km
  T = hdf5read(InFile, 't_coords') / 3600; % hr

  % Reduce data to single Z level
  Pdata = squeeze(Hvar(:,:,Zlevel,:));
  clear Hvar;
  
  % Variable will be of the form: (Lon, Lat, Time)
  [ Nlon, Nlat, Nt ] = size(Pdata);
 
  % Create the movie - go frame by frame creating a plot each time and
  % writing that plot into the output file.
  fprintf('Writing GIF movie file: %s\n', OutFile);

  Fig = figure;
  
  % get the frame size
  Pos = get(gcf, 'Position');
  Width = Pos(3);
  Height = Pos(4);
  
  % preallocate the movie buffer
  Tsteps = Tstart:TsInc:Tend;
  mov = zeros(Height, Width, 1, length(Tsteps), 'uint8');

  for i = Tsteps
    fprintf('  Creating frame for time step: %d\n', i);
    
    % make the plot
    contourf(Lon,Lat,Pdata(:,:,i)');
    caxis([ Pmin Pmax ]);
    daspect(Paspect);
    shading flat;
    colorbar;
    colormap('cool');
    title(Ptitle);
    
    % create the frame
    drawnow;
    frame = getframe(1);
    im = frame2im(frame);
    if (i == 1)
      [ mov(:,:,1,i), cmap ] = rgb2ind(im, 256);
    else
      mov(:,:,1,i) = rgb2ind(im, cmap);
    end
  end

  close(Fig);
  
  imwrite(mov, cmap, OutFile,'gif', 'DelayTime', Delay, 'Loopcount',inf);
  
end
