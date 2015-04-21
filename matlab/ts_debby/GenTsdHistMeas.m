function [ ] = GenTsdHistMeas(ConfigFile)

  [ Config ] = ReadConfig(ConfigFile);

  Adir = Config.AzavgDir;
  Ddir = Config.DiagDir;

  % make sure output directory exists
  if (exist(Ddir, 'dir') ~= 7)
    mkdir(Ddir);
  end

  % Description of measurements
  %   <measure_name> <measure_list> <out_file_prefix>
  %
  %   where <measure_list> is one or more of:
  %     <file_prefix> <var_name> <reduction_method> <arg_for_reduction_method> <out_var_name> <select_op> <select_val>
  MeasSets = {
    % storm speed measurements
%    {
%      'Speed'
%      {
%        { 'hist_speed' 'speed' 'farea' 0.50 'avg_speed'     'ge' 0   }
%        { 'hist_speed' 'speed' 'farea' 0.95 'max_speed'     'ge' 0   }
%  
%        { 'hist_speed10m' 'speed10m' 'farea' 0.50 'avg_speed10m'     'ge' 0   }
%        { 'hist_speed10m' 'speed10m' 'farea' 0.95 'max_speed10m'     'ge' 0   }
%  
%        { 'hist_speed_t' 'speed_t' 'farea' 0.50 'avg_speed_t'     'ge' 0   }
%        { 'hist_speed_t' 'speed_t' 'farea' 0.95 'max_speed_t'     'ge' 0   }
%  
%        { 'hist_speed_r' 'speed_r' 'farea' 0.50 'avg_speed_r'     'ge' 0   }
%        { 'hist_speed_r' 'speed_r' 'farea' 0.95 'max_speed_r'     'ge' 0   }
%      }
%      'hist_meas_speed'
%    }
%  
%    % storm pressure measurements
%    {
%      'Pressure'
%      {
%        { 'hist_press' 'press' 'farea' 0.50 'avg_press'     'ge' 0   }
%        { 'hist_press' 'press' 'farea' 0.05 'min_press'     'ge' 0   }
%  
%        { 'hist_sea_press' 'sea_press' 'farea' 0.50 'avg_sea_press'     'ge' 0   }
%        { 'hist_sea_press' 'sea_press' 'farea' 0.05 'min_sea_press'     'ge' 0   }
%      }
%      'hist_meas_press'
%    }
%
%    % precip rate measurements
%    {
%      'Precip Rate'
%      {
%        { 'hist_pcprate' 'pcprate' 'farea' 0.50 'avg_pcprate'     'ge' 0   }
%        { 'hist_pcprate' 'pcprate' 'farea' 0.95 'max_pcprate'     'ge' 0   }
%      }
%      'hist_meas_pcprate'
%    }

    % vertical velocity measurements
    {
      'Vertical Velocity'
      {
        { 'hist_w' 'w' 'farea' 0.50 'avg_updraft'     'ge'  0.1   }
        { 'hist_w' 'w' 'farea' 0.95 'max_updraft'     'ge'  0.1   }

        { 'hist_w' 'w' 'farea' 0.50 'avg_dndraft'     'le' -0.1   }
        { 'hist_w' 'w' 'farea' 0.05 'max_dndraft'     'le' -0.1   }
      }
      'hist_meas_w'
    }


    };

  Nsets = length(MeasSets);

  for icase = 1:length(Config.Cases)
    Case = Config.Cases(icase).Cname;

    fprintf('*****************************************************************\n');
    fprintf('Generating histogram measurements:\n');
    fprintf('  Case: %s\n', Case);
    fprintf('\n');


    for iset = 1:Nsets
      MeasName   = MeasSets{iset}{1};
      MeasList   = MeasSets{iset}{2};
      OutFprefix = MeasSets{iset}{3};

      fprintf('    Measurement set: %s\n', MeasName);
      fprintf('\n');

      % Put all measurements into one file per case
      % If the file exists, remove it so that the HDF5 commands
      % can effectively re-create datasets.
      OutFile = sprintf('%s/%s_%s.h5', Ddir, OutFprefix, Case);
      if (exist(OutFile, 'file') == 2)
        delete(OutFile);
      end
  
      Nmeas = length(MeasList);
      MaxNz = 0;
      for imeas = 1:Nmeas
        Fprefix   = MeasList{imeas}{1};
        Vname     = MeasList{imeas}{2};
        Rmethod   = MeasList{imeas}{3};
        Param     = MeasList{imeas}{4};
        OutVname  = MeasList{imeas}{5};
        SelectOp  = MeasList{imeas}{6};
        SelectVal = MeasList{imeas}{7};
  
        % add on leading '/' for HDF5 routines
        Vname    = sprintf('/%s', Vname);
        OutVname = sprintf('/%s', OutVname);
  
        InFile = sprintf('%s/%s_%s.h5', Adir, Fprefix, Case);
  
        fprintf('      Reading: %s (%s)\n', InFile, Vname);
        fprintf('        Reduction method: %s (%.2f)\n', Rmethod, Param);
        fprintf('        Selection: %s %.2f\n', SelectOp, SelectVal);
  
        % Read in data which will be 4D -> (x,y,z,t)
        %
        %     x --> radial bands
        %     y --> histogram bins
        %     z --> height
        %     t --> time
        %
        HDATA = squeeze(h5read(InFile, Vname));
        BINS  = squeeze(h5read(InFile, '/y_coords'));
  
        % Assume same r,z,t values for all measurements
        X    = squeeze(h5read(InFile, '/x_coords'));
        Y    = 1; % dummy dimension for output
        Z    = squeeze(h5read(InFile, '/z_coords'));
        T    = squeeze(h5read(InFile, '/t_coords'));
  
        Nx = length(X);
        Ny = length(Y);
        Nz = length(Z);
        Nt = length(T);

        % Keep track of the coordinates of the full (maximum sized) z dimension
        if (Nz > MaxNz)
          MaxNz = Nz;
          AllNz = Nz;
          AllZ = Z;
        end
  
        % Reduce the histograms level by level to preserve vertical structure
        % MEAS will be organized as: (r,z,t)
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
  
        MEAS = squeeze(ReduceHists(HDATA(:,B1:B2,:,:), 2, BINS(B1:B2), Rmethod, Param));
  
        % Write out measurement
        fprintf('      Writing: %s (%s)\n', OutFile, OutVname)
        fprintf('\n');
  
        % Write out measurement -> force to be (x,y,z,t) for dimension
        % attach code below.
        OutVar = reshape(MEAS, [ Nx Ny Nz Nt ]);
        h5create(OutFile, OutVname, size(OutVar));
        h5write(OutFile, OutVname, OutVar);
      end % measurements
  
      % Create the dimensions
      fprintf('      Creating dimensions: %s\n', OutFile);
      fprintf('\n');
  
      Xname = '/x_coords';
      Yname = '/y_coords';
      Zname = '/z_coords';
      Tname = '/t_coords';
  
      CreateDimensionsXyzt(OutFile, X, Y, AllZ, T, Xname, Yname, Zname, Tname);
  
      % Attach dimensions to all variables
      for imeas = 1:Nmeas
        Vname = MeasList{imeas}{5};   % use the output var name
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
    end % sets
  end % cases
end % function
