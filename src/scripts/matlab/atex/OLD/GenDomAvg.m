% Script to generate domain averaged data.
%
% This script assumes that the time steps go in order when you
% concatenate together the results from each time directory
% in the order specified in Tdirs.
%
% This allows this script to build an array with each profile
% in a row and then call the mean() function to generated the
% overall avg profile.

clear;

VarLH = 'lh_vapt';
HvarLH = '/latheatvapt';
VarSW = 'swup';
HvarSW = '/swup';
OutDir = 'DIAGS';

% Create the output directory if it doesn't exist
if (exist(OutDir,'dir') == 0)
  mkdir(OutDir);
end

[ Exps, Tdirs ] = ReadConfig('Config');

fprintf('Running domain averaging on: %s\n', VarLH);
for i_exp = 1:length(Exps)
  fprintf('  Experiment: %s\n', char(Exps(i_exp)));

  for i_tdir = 1:length(Tdirs)
    % Form the directory: <Exp>/HDF5/<Tdir>
    LHdir = sprintf('%s/HDF5/%s', char(Exps(i_exp)), char(Tdirs(i_tdir)));
    LHfpat = sprintf('%s/%s-*', LHdir, VarLH);
    LHfile = dir(LHfpat);
    LHfile = LHfile.name;

    SWdir = sprintf('%s/HDF5/%s', char(Exps(i_exp)), char(Tdirs(i_tdir)));
    SWfpat = sprintf('%s/%s-*', SWdir, VarSW);
    SWfile = dir(SWfpat);
    SWfile = SWfile.name;
    
    % Need to pre-pend the directory since dir returns the basename of the file
    LHfile = sprintf('%s/%s', LHdir, LHfile);
    fprintf('    Reading HDF5 file: %s\n', LHfile);

    % Read in data for the variable from this time directory
    LHdata = hdf5read(LHfile, HvarLH);
    Zcoords = hdf5read(LHfile, '/z_coords');
    Tcoords = hdf5read(LHfile, '/t_coords');

    % Need to pre-pend the directory since dir returns the basename of the file
    SWfile = sprintf('%s/%s', SWdir, SWfile);
    fprintf('    Reading HDF5 file: %s\n', SWfile);

    % Read in data for the variable from this time directory
    SWdata = hdf5read(SWfile, HvarSW);

    % Trim off the lateral boundaries since these are set to zeros
    %  Vars are: (x,y,z,t)
    LHdata = LHdata(2:end-1,2:end-1,:,:);
    SWdata = SWdata(2:end-1,2:end-1,:,:);

    % Find the average for this piece
    % Vars are: (x,y,z,t)
    DomAvgLH  = squeeze(mean(mean(LHdata,1),2));
    DomAvgSW  = squeeze(mean(mean(SWdata,1),2));

    % Put this into a 2D array: (z,t)
    if (i_tdir == 1)
      AllDomAvgLH = DomAvgLH;
      AllDomAvgSW = DomAvgSW;
      AllTcoords = Tcoords;
    else
      AllDomAvgLH = cat(2, AllDomAvgLH, DomAvgLH);
      AllDomAvgSW = cat(2, AllDomAvgSW, DomAvgSW);
      AllTcoords = cat(1, AllTcoords, Tcoords);
    end
  end
  fprintf('\n');

  % Save the averaged data in an hdf5 file
  LHfile = sprintf('./%s/%s_tdavg_%s.h5', OutDir, VarLH, char(Exps(i_exp)));
  Hdset = sprintf('/%s',VarLH);
  fprintf('    Writing HDF5 file: %s\n', LHfile);
  hdf5write(LHfile, Hdset, AllDomAvgLH);
  hdf5write(LHfile, '/z_coords', Zcoords, 'WriteMode', 'append');
  hdf5write(LHfile, '/t_coords', AllTcoords, 'WriteMode', 'append');
  fprintf('\n');

  SWfile = sprintf('./%s/%s_tdavg_%s.h5', OutDir, VarSW, char(Exps(i_exp)));
  Hdset = sprintf('/%s',VarSW);
  fprintf('    Writing HDF5 file: %s\n', SWfile);
  hdf5write(SWfile, Hdset, AllDomAvgSW);
  hdf5write(SWfile, '/z_coords', Zcoords, 'WriteMode', 'append');
  hdf5write(SWfile, '/t_coords', AllTcoords, 'WriteMode', 'append');
  fprintf('\n');
end
