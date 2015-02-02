function [ ] = GenDiffProfiles(ConfigFile)

  [ Config ] = ReadConfig(ConfigFile);

  Ddir = Config.DiagDir;

  % List of cases, put control case first
  CaseList = {
    'TSD_MOIST_NODUST'
    'TSD_MOIST_DUST'
    'TSD_DRY_NODUST'
    'TSD_DRY_DUST'
    };
  Ncase = length(CaseList);

  % Description of profiles
  %   <file_prefix> <var_name> <reduction_method> <arg_for_reduction_method> <out_var_name> <select_op> <select_val>
  ProfList = {
    { 'hist_profs' 'prof_ccn_conc_tavg' 'ccn_conc_tavg' 0 }
    { 'hist_profs' 'prof_d1_conc_tavg'  'd1_conc_tavg'  0 }
    { 'hist_profs' 'prof_d2_conc_tavg'  'd2_conc_tavg'  0 }

    { 'hist_profs' 'prof_lat_heat_vapt_tavg' 'lat_heat_vapt_tavg' 0 }
    { 'hist_profs' 'prof_lat_cool_vapt_tavg' 'lat_cool_vapt_tavg' 0 }

    { 'hist_profs' 'prof_lat_heat_frzt_tavg' 'lat_heat_frzt_tavg' 0 }
    { 'hist_profs' 'prof_lat_cool_frzt_tavg' 'lat_cool_frzt_tavg' 0 }

    { 'hist_profs' 'prof_w_up_tavg'   'w_up_tavg'   0 }
    { 'hist_profs' 'prof_w_down_tavg' 'w_down_tavg' 0 }

    { 'hist_profs' 'prof_theta_e_tavg' 'theta_e_tavg' 0 }
    };
    Nprof = length(ProfList);

  % Put all profiles into one file
  % If the file exists, remove it so that the HDF5 commands
  % can effectively re-create datasets.
  OutFile = sprintf('%s/diff_profs.h5', Ddir);
  if (exist(OutFile, 'file') == 2)
    delete(OutFile);
  end

  fprintf('*****************************************************************\n');
  fprintf('Generating difference profiles:\n');
  fprintf('  Output File: %s\n', OutFile);
  fprintf('\n');


  % For each profile, walk through cases. Record first case as control, then for
  % other cases, subtract off the control profile.
  for iprof = 1:Nprof
    Fprefix   = ProfList{iprof}{1};
    Vname     = ProfList{iprof}{2};
    OutVname  = ProfList{iprof}{3};
    NanVal    = ProfList{iprof}{4};

    fprintf('  Profile: %s\n', Vname);

    % add on leading '/' for HDF5 routines
    Vname    = sprintf('/%s', Vname);
    OutVname = sprintf('/%s', OutVname);

    for icase = 1:Ncase
      Case = CaseList{icase};
      fprintf('  Case: %s\n', Case);

      InFile = sprintf('%s/%s_%s.h5', Ddir, Fprefix, Case);

      if (icase == 1)
        fprintf('    Reading (control): %s (%s)\n', InFile, Vname);
        CONTROL_PROF = squeeze(h5read(InFile, Vname));
        CONTROL_PROF(isnan(CONTROL_PROF)) = NanVal;

        X    = squeeze(h5read(InFile, '/x_coords'));
        Y    = squeeze(h5read(InFile, '/y_coords'));
        Z    = squeeze(h5read(InFile, '/z_coords'));
        T    = squeeze(h5read(InFile, '/t_coords'));

        Nx = length(X);
        Ny = length(Y);
        Nz = length(Z);
        Nt = length(T);
      else
        fprintf('    Reading: %s (%s)\n', InFile, Vname);
        PROF = squeeze(h5read(InFile, Vname));
        PROF(isnan(PROF)) = NanVal;

        DIFF_PROF = PROF - CONTROL_PROF;

        % Write out difference profile
        OutDset = sprintf('%s_%s', OutVname, Case);
        h5create(OutFile, OutDset, size(DIFF_PROF));
        h5write(OutFile, OutDset, DIFF_PROF);
      end

    end % cases
    fprintf('\n');
  end   % profiles

  % Write out the coordinate values
  h5create(OutFile, '/x_coords', Nx);
  h5write(OutFile, '/x_coords', X);

  h5create(OutFile, '/y_coords', Ny);
  h5write(OutFile, '/y_coords', Y);

  h5create(OutFile, '/z_coords', Nz);
  h5write(OutFile, '/z_coords', Z);

  h5create(OutFile, '/t_coords', Nt);
  h5write(OutFile, '/t_coords', T);

  % Now spread out the difference profiles into separate files (case by case) so that
  % GenContourPlots() will be able to find them.
  InFile = sprintf('%s/diff_profs.h5', Ddir);
  fprintf('Creating case specific output files:\n');
  for icase = 2:Ncase
    Case = CaseList{icase};
    fprintf('  Case: %s\n', Case);

    OutFile = sprintf('%s/diff_profs_%s.h5', Ddir, Case);
    if (exist(OutFile, 'file') == 2)
      delete(OutFile);
    end

    for iprof = 1:Nprof
      Vname = ProfList{iprof}{3};

      InVname = sprintf('/%s_%s', Vname, Case);
      OutVname = sprintf('/%s', Vname);

      fprintf('    Reading: %s (%s)\n', InFile, InVname);

      PROF = squeeze(h5read(InFile, InVname));

      fprintf('    Writing: %s (%s)\n', OutFile, OutVname);
      h5create(OutFile, OutVname, size(PROF));
      h5write(OutFile, OutVname, PROF);
    end

    % Copy over the coordinate values
    X    = squeeze(h5read(InFile, '/x_coords'));
    Y    = squeeze(h5read(InFile, '/y_coords'));
    Z    = squeeze(h5read(InFile, '/z_coords'));
    T    = squeeze(h5read(InFile, '/t_coords'));

    h5create(OutFile, '/x_coords', Nx);
    h5write(OutFile, '/x_coords', X);

    h5create(OutFile, '/y_coords', Ny);
    h5write(OutFile, '/y_coords', Y);

    h5create(OutFile, '/z_coords', Nz);
    h5write(OutFile, '/z_coords', Z);

    h5create(OutFile, '/t_coords', Nt);
    h5write(OutFile, '/t_coords', T);

    fprintf('\n');
  end
end     % function
