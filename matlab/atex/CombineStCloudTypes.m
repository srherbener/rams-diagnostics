function [ ] = CombineStCloudTypes(ConfigFile)
% CombineStCloudTypes combine sum terms and npts for strat, strnp, stmix data

  % Read the config file to get the structure of how the data is laid out in
  % the file system.
  [ Config ] = ReadConfig(ConfigFile);

  Ddir = Config.DiagDir;

  % Format for FileSets entry
  %   { <list_of_input_files> <list_of_variables> <reshape_to_4d> <output_file> }
  FileSets = {
    { { 'w_theta_flux_strat'   'w_theta_flux_stmix'   'w_theta_flux_strnp'   } { 'num_points' 'w-theta'   } 'w_theta_flux_stall'    }
    { { 'w_theta_v_flux_strat' 'w_theta_v_flux_stmix' 'w_theta_v_flux_strnp' } { 'num_points' 'w-theta_v' } 'w_theta_v_flux_stall'  }
    { { 'w_speed_flux_strat'   'w_speed_flux_stmix'   'w_speed_flux_strnp'   } { 'num_points' 'w-speed'   } 'w_speed_flux_stall'    }
    { { 'w_vapor_flux_strat'   'w_vapor_flux_stmix'   'w_vapor_flux_strnp'   } { 'num_points' 'w-vapor'   } 'w_vapor_flux_stall'    }
    };

  for icase = 1:length(Config.Cases)
    Case = Config.Cases(icase).Cname;
    fprintf('**************************************************************\n');
    fprintf('Combining St data for case: %s\n', Case);
    fprintf('\n');
 
    % for each file set, walk through each input file and accumulate
    for ifset = 1:length(FileSets)
      InFprefixList = FileSets{ifset}{1};
      InVarList     = FileSets{ifset}{2};
      OutFprefix    = FileSets{ifset}{3};

      fprintf('  Variables:\n');
      Nvar = length(InVarList);
      for ivar = 1:Nvar
        fprintf('    %s\n', InVarList{ivar});
      end
      fprintf('\n');

      fprintf('  Files\n');
      Nfile = length(InFprefixList);
      for ifile = 1:Nfile
        InFile = sprintf('%s/%s_%s.h5', Ddir, InFprefixList{ifile}, Case);
        fprintf('    %s\n', InFile);

        for ivar = 1:Nvar
          HDATA = hdf5read(InFile, InVarList{ivar});

          if (ifile == 1)
            VAR_SUMS{ivar} = HDATA;

            if (ivar == 1)
              X = h5read(InFile, '/x_coords');
              Y = h5read(InFile, '/y_coords');
              Z = h5read(InFile, '/z_coords');
              T = h5read(InFile, '/t_coords');
            end
          else
             VAR_SUMS{ivar} = VAR_SUMS{ivar} + HDATA;
          end
        end
      end
      fprintf('\n');

      OutFile = sprintf('%s/%s_%s.h5', Ddir, OutFprefix, Case);
      fprintf('  Writing: %s\n', OutFile);

      hdf5write(OutFile, '/header', 'ATEX');
      for ivar = 1:Nvar
        hdf5write(OutFile, InVarList{ivar}, VAR_SUMS{ivar}, 'WriteMode', 'append');
      end

      hdf5write(OutFile, '/x_coords', X, 'WriteMode', 'append');
      hdf5write(OutFile, '/y_coords', Y, 'WriteMode', 'append');
      hdf5write(OutFile, '/z_coords', Z, 'WriteMode', 'append');
      hdf5write(OutFile, '/t_coords', T, 'WriteMode', 'append');
      fprintf('\n');
    end
  end
end
