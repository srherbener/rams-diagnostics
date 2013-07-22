function [ ] = CompileCounts(ConfigFile)
% CompileCounts compile count data for subsequent POP slope calculation

% Read the config file to get the structure of how the data is laid out in
% the file system.
[ Config ] = ReadConfig(ConfigFile);

Tdir = Config.TsavgDir;
Ddir = Config.DiagDir;

% Make sure output directory exists
if (exist(Ddir, 'dir') ~= 7)
    mkdir(Ddir);
end

for ismeas = 1: length(Config.CompCounts)
  % clear these out because they may need to change size, and without the clear
  % they will just retain the largest size (and associated data) that they have
  % been set to so far
  clear Xvals;
  clear ALL_COUNTS;;

  Name = Config.CompCounts(ismeas).Name;
  InDir = Config.CompCounts(ismeas).InDir;
  Fprefix = Config.CompCounts(ismeas).Fprefix;
  Vname = Config.CompCounts(ismeas).Rvar;

  fprintf('***********************************************************************\n');
  fprintf('Generating Slope time series:\n');
  fprintf('  Name: %s\n', Name);
  fprintf('  Variable: %s\n', Vname);
  fprintf('\n');

  OutFile = sprintf('%s/%s.h5', Ddir, Name);

  ips = Config.CompCounts(ismeas).PSnum;
  if (ips == 0)
    fprintf('  WARNING: skipping CompCounts number %d due to no associated PlotSet\n', ismeas)
  else
    % For each case in the plot set, read in the count data and simply add together
    % the arrays of counts from the group of files.
    for icase = 1:Config.PlotSets(ips).Ncases
      Case = Config.PlotSets(ips).Cases(icase).Cname;

      InFile  = sprintf('%s/%s_%s.h5', InDir, Fprefix, Case);
      Hdset   = sprintf('/%s', Vname);

      % Read in the count data. COUNTS will be organized as (x,y,z,t) where
      %    x --> LWP bins
      %    y --> LTSS bins
      %    z --> counts: 1 -> Nt, 2 -> Nr
      %    t --> time
      fprintf('  Reading file: %s, Dataset: %s\n', InFile, Hdset);
      COUNTS = hdf5read(InFile, Hdset);
      if (icase == 1)
        % if on the first file, read in the coordinates and load up the
        % output array with COUNTS.
        X = hdf5read(InFile, '/x_coords');
        Y = hdf5read(InFile, '/y_coords');
        Z = hdf5read(InFile, '/z_coords');
        T = hdf5read(InFile, '/t_coords');

        ALL_COUNTS = COUNTS;
      else
        % on 2..n files --> add in COUNTS to ALL_COUNTS.
        ALL_COUNTS = ALL_COUNTS + COUNTS;
      end
    end
    fprintf('\n');

    fprintf('  Writing file: %s\n', OutFile);
    Hdset = sprintf('/%s', Vname);
    hdf5write(OutFile, Hdset, ALL_COUNTS);

    hdf5write(OutFile, '/x_coords', X, 'WriteMode', 'append');
    hdf5write(OutFile, '/y_coords', Y, 'WriteMode', 'append');
    hdf5write(OutFile, '/z_coords', Z, 'WriteMode', 'append');
    hdf5write(OutFile, '/t_coords', T, 'WriteMode', 'append');

    % attach the dimensions
    fprintf('    Attaching dimensions: %s\n', OutFile);
    file_id = H5F.open(OutFile, 'H5F_ACC_RDWR', 'H5P_DEFAULT');
    var_dset_id = H5D.open(file_id, Hdset, 'H5P_DEFAULT');
    x_dset_id = H5D.open(file_id, 'x_coords', 'H5P_DEFAULT');
    y_dset_id = H5D.open(file_id, 'y_coords', 'H5P_DEFAULT');
    z_dset_id = H5D.open(file_id, 'z_coords', 'H5P_DEFAULT');
    t_dset_id = H5D.open(file_id, 't_coords', 'H5P_DEFAULT');
  
    % set x,y,z,t vars as dimensions
    H5DS.set_scale(x_dset_id, 'x');
    H5DS.set_scale(y_dset_id, 'y');
    H5DS.set_scale(z_dset_id, 'z');
    H5DS.set_scale(t_dset_id, 't');
  
    % attach dimensions to
    % the dimension order gets reversed (third argument to
    % the H5DS.attach_scale routine) since the dimensions 
    % get reversed in the HDF5 file.
    H5DS.attach_scale(var_dset_id, x_dset_id, 3);
    H5DS.attach_scale(var_dset_id, y_dset_id, 2);
    H5DS.attach_scale(var_dset_id, z_dset_id, 1);
    H5DS.attach_scale(var_dset_id, t_dset_id, 0);
  
    % clean up
    H5D.close(x_dset_id);
    H5D.close(y_dset_id);
    H5D.close(z_dset_id);
    H5D.close(t_dset_id);
    H5D.close(var_dset_id);
    H5F.close(file_id);

    fprintf('\n');
  end
end

end
