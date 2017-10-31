function [] = GenRadarData(ConfigFile)
% GenRadarData fucntion to generate radar data (horizontal slice from 3D reflectivity field)

[ Config ] = ReadConfig(ConfigFile);

Ddir = Config.DiagDir;

% Make sure output directory exists
if (exist(Ddir, 'dir') ~= 7)
    mkdir(Ddir);
end

Xmin = -180;
Xmax = 180;
Ymin = -180;
Ymax = 180;
Zmin = 1800;     % m AGL
Zmax = 2000;     % horizontal slice at 1970 (near 2000)
Tmin = 24*3600;    % sec
Tmax = 144*3600;   % sec

Hvar = 'reflect_all';

for icase = 1:length(Config.Cases)
  Case = Config.Cases(icase).Cname;

  Hfile = sprintf('HDF5/%s-%s-AS-1998-08-22-120000-g3.h5', Hvar, Case);
  OutFile = sprintf('%s/%s_%s.h5', Ddir, Hvar, Case);

  fprintf('***********************************************************************\n');
  fprintf('Generating radar reflectivity data:\n');
  fprintf('  Case: %s\n', Case);
  fprintf('  Variable: %s\n', Hvar);
  fprintf('  Input file: %s\n', Hfile);
  fprintf('  Output file: %s\n', OutFile);
  fprintf('  Sample level: %d m\n', Zmin);
  fprintf('\n');

  X = ReadSelectXyzt(Hfile, 'x_coords', Xmin, Xmax, Ymin, Ymax, Zmin, Zmax, Tmin, Tmax);
  Y = ReadSelectXyzt(Hfile, 'y_coords', Xmin, Xmax, Ymin, Ymax, Zmin, Zmax, Tmin, Tmax);
  T = ReadSelectXyzt(Hfile, 't_coords', Xmin, Xmax, Ymin, Ymax, Zmin, Zmax, Tmin, Tmax);
  RDBZ = ReadSelectXyzt(Hfile, Hvar, Xmin, Xmax, Ymin, Ymax, Zmin, Zmax, Tmin, Tmax);
  RDBZ(RDBZ < -2.5) = -2.5; % wipe out dBZ less than -2.5 so the plotting colors will work out like NWS radar images

  Nx = length(X);
  Ny = length(Y);
  Nt = length(T);

  % RDBZ is now (x,y,t), need to insert a dummy z dimension
  Z = 0;
  OutVar = reshape(RDBZ, [ Nx Ny 1 Nt ]);

  hdf5write(OutFile, Hvar, OutVar);
  hdf5write(OutFile, '/x_coords', X, 'WriteMode', 'append');
  hdf5write(OutFile, '/y_coords', Y, 'WriteMode', 'append');
  hdf5write(OutFile, '/z_coords', Z, 'WriteMode', 'append');
  hdf5write(OutFile, '/t_coords', T, 'WriteMode', 'append');
  
  fprintf('\n');
end

end
