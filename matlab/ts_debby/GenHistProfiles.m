function [ ] = GenHistProfiles(ConfigFile)

  [ Config ] = ReadConfig(ConfigFile);

  Adir = Config.AzavgDir;
  Ddir = Config.DiagDir;

  % make sure output directory exists
  if (exist(Ddir, 'dir') ~= 7)
    mkdir(Ddir);
  end

  % Description of profiles
  %   <file_prefix> <var_name> <reduction_method> <arg_for_reduction_method> <out_var_name> <select_op> <select_val>
  ProfList = {
    { 'hist_ccn_conc' 'ccn_conc' 'wtmean' 0.5 'prof_ccn_conc'      'ge' 20   }
    { 'hist_d1_conc'  'd1_conc'  'wtmean' 0.5 'prof_d1_conc'       'ge' 20   }
    { 'hist_d2_conc'  'd2_conc'  'wtmean' 0.5 'prof_d2_conc'       'ge' 20   }

    { 'hist_lh_vapt'  'lh_vapt'  'wtmean' 0.5 'prof_lat_heat_vapt' 'ge'  10  }
    { 'hist_lh_vapt'  'lh_vapt'  'wtmean' 0.5 'prof_lat_cool_vapt' 'le' -10  }

    { 'hist_lh_frzt'  'lh_frzt'  'wtmean' 0.5 'prof_lat_heat_frzt' 'ge'  1   }
    { 'hist_lh_frzt'  'lh_frzt'  'wtmean' 0.5 'prof_lat_cool_frzt' 'le' -1   }

    { 'hist_w'        'w'        'wtmean' 0.5 'prof_w_up'          'ge'  0.1 }
    { 'hist_w'        'w'        'wtmean' 0.5 'prof_w_down'        'le' -0.1 }

    { 'hist_theta_e'  'theta_e'  'wtmean' 0.5 'prof_theta_e'       'ge'  0   }
    };

    Nprof = length(ProfList);

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

    for iprof = 1:Nprof
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
        X    = squeeze(h5read(InFile, '/x_coords'));
        Y    = 1; % dummy dimension for output
        Z    = squeeze(h5read(InFile, '/z_coords'));
        T    = squeeze(h5read(InFile, '/t_coords'));

        SIM_TIME = (T ./ 3600) - 42; % h, starting with 0
        T1 = find(SIM_TIME >= 10, 1, 'first');

        Nx = length(X);
        Ny = length(Y);
        Nz = length(Z);
        Nt = length(T);
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
      TAVG_PROF = nanmean(PROF(:,:,T1:end),3);

      % Write out measurement
      fprintf('  Writing: %s (%s)\n', OutFile, OutVname)

      % Write out profile -> force to be (x,y,z,t) for dimension
      % attach code below.
      OutVar = reshape(PROF, [ Nx Ny Nz Nt ]);
      h5create(OutFile, OutVname, size(OutVar));
      h5write(OutFile, OutVname, OutVar);

      OutVar = reshape(TAVG_PROF, [ Nx Ny Nz 1 ]);
      TavgOutVname = sprintf('%s_tavg', OutVname);
      h5create(OutFile, TavgOutVname, size(OutVar));
      h5write(OutFile, TavgOutVname, OutVar);
    end

    % Create the dimensions
    fprintf('  Creating dimensions: %s\n', OutFile);
    Xname = '/x_coords';
    Yname = '/y_coords';
    Zname = '/z_coords';
    Tname = '/t_coords';

    CreateDimensionsXyzt(OutFile, X, Y, Z, T, Xname, Yname, Zname, Tname);

    % Attach dimensions to all variables
    for iprof = 1:Nprof
      Vname = ProfList{iprof}{5};   % use the output var name
      AttachDimensionsXyzt(OutFile, Vname, Xname, Yname, Zname, Tname);
    end

    % GRADS needs the following attributes on the dimension datasets in order
    % to recognize which dimensions go with which datasets. These attribute
    % names and values are following the COARDS convention.
    WriteStringAttribute(OutFile, Xname, 'axis', 'x');
    WriteStringAttribute(OutFile, Xname, 'units', 'degrees_east');

    WriteStringAttribute(OutFile, Yname, 'axis', 'y');
    WriteStringAttribute(OutFile, Yname, 'units', 'degrees_north');

    WriteStringAttribute(OutFile, Zname, 'axis', 'z');
    WriteStringAttribute(OutFile, Zname, 'units', 'meters');

    WriteStringAttribute(OutFile, Tname, 'axis', 't');
    WriteStringAttribute(OutFile, Tname, 'units', 'seconds since 2006-08-20 12:00:00 00:00');
    
    fprintf('\n');
  end
end 
