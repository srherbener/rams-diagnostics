function [ ] = GenHistProfiles(ConfigFile)

  [ Config ] = ReadConfig(ConfigFile);

  Adir = Config.AzavgDir;
  Ddir = Config.DiagDir;

  % make sure output directory exists
  if (exist(Ddir, 'dir') ~= 7)
    mkdir(Ddir);
  end

  % Description of profiles
  ProfList = {
    { 'hist_ccn_conc' 'ccn_conc' 'wtmean' 0.5 'prof_ccn_conc' 'ge' 20 }
    { 'hist_d1_conc'  'd1_conc'  'wtmean' 0.5 'prof_d1_conc'  'ge' 20 }
    { 'hist_d2_conc'  'd2_conc'  'wtmean' 0.5 'prof_d2_conc'  'ge' 20 }

    { 'hist_lh_vapt'  'lh_vapt'  'wtmean' 0.5 'prof_lat_heat_vapt' 'ge'  10 }
    { 'hist_lh_vapt'  'lh_vapt'  'wtmean' 0.5 'prof_lat_cool_vapt' 'le' -10 }

    { 'hist_lh_frzt'  'lh_frzt'  'wtmean' 0.5 'prof_lat_heat_frzt' 'ge'  1 }
    { 'hist_lh_frzt'  'lh_frzt'  'wtmean' 0.5 'prof_lat_cool_frzt' 'le' -1 }
    };

  for icase = 1:length(Config.Cases)
    Case = Config.Cases(icase).Cname;

    fprintf('*****************************************************************\n');
    fprintf('Generating profiles:\n');
    fprintf('  Case: %s\n', Case);
    fprintf('\n');

    % Put all profiles into one file per case
    % If the file exists, remove it so that the HDF5 commands
    % can effectively re-create datasets.
    OutFile = sprintf('%s/hist_profs_%s.h5', Ddir, Case);
    if (exist(OutFile, 'file') == 2)
      delete(OutFile);
    end

    for iprof = 1:length(ProfList)
      Fprefix   = ProfList{iprof}{1};
      Vname     = ProfList{iprof}{2};
      Rmethod   = ProfList{iprof}{3};
      Ptile     = ProfList{iprof}{4};
      OutVname  = ProfList{iprof}{5};
      SelectOp  = ProfList{iprof}{6};
      SelectVal = ProfList{iprof}{7};

      % add on leading '/' for HDF5 routines
      Vname    = sprintf('/%s', Vname);
      OutVname = sprintf('/%s', OutVname);

      InFile = sprintf('%s/%s_%s.h5', Adir, Fprefix, Case);

      fprintf('  Reading: %s (%s)\n', InFile, Vname);
      fprintf('    Reduction method: %s (%.2f)\n', Rmethod, Ptile);
      fprintf('    Selection: %s %.2f\n', SelectOp, SelectVal);

      % Read in data which will be 4D -> (x,y,z,t)
      %
      %     x --> radial bands
      %     y --> histogram bins
      %     z --> height
      %     t --> time
      %
      HDATA = squeeze(h5read(InFile, Vname));
      BINS  = squeeze(h5read(InFile, '/y_coords'));

      % Assume same r,z,t values for all profiles
      if (iprof == 1)
        R    = squeeze(h5read(InFile, '/x_coords'));
        Z    = squeeze(h5read(InFile, '/z_coords'));
        T    = squeeze(h5read(InFile, '/t_coords'));
      end

      % Reduce the histograms level by level to create vertical profiles
      % PROF will be organized as: (r,z,t)
      %
      % Do selection on bin values (ge or le)
      %

      % select all bin values by default
      B1 = 1;         
      B2 = length(BINS);

      if (strcmp(SelectOp, 'ge'))
        B1 = find(BINS >= SelectVal, 1, 'first');
      end

      if (strcmp(SelectOp, 'le'))
        B2 = find(BINS <= SelectVal, 1, 'last');
      end

      PROF = squeeze(ReduceHists(HDATA(:,B1:B2,:,:), 2, BINS(B1:B2), Rmethod, Ptile));

      % Write out measurement
      fprintf('  Writing: %s (%s)\n', OutFile, OutVname)

      % Write out profile
      h5create(OutFile, OutVname, size(PROF));
      h5write(OutFile, OutVname, PROF);
    end

    % write out the coordinate values
    h5create(OutFile, '/x_coords', size(R));
    h5write(OutFile, '/x_coords', R);

    h5create(OutFile, '/z_coords', size(Z));
    h5write(OutFile, '/z_coords', Z);

    h5create(OutFile, '/t_coords', size(T));
    h5write(OutFile, '/t_coords', T);

    fprintf('\n');
  end
end 
