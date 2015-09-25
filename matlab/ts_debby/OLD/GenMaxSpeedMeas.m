function [ ] = GenMaxSpeedMeas()

  Ddir = 'DIAGS';
  if (exist(Ddir, 'dir') ~= 7)
    mkdir(Ddir);
  end

  Cases = {
   'TSD_SAL_DUST'
   'TSD_SAL_NODUST'
   'TSD_NONSAL_DUST'
   'TSD_NONSAL_NODUST'
   };
  Ncases = length(Cases);

  ControlCase = 'TSD_SAL_DUST';

  % Description of measurements
  MeasList = {
    % in_file in_var out_var sfc_out_var
    { 'DIAGS/hist_meas_speed_<CASE>.h5' '/ps_max_speed'   '/ps_max_speed'   '/ps_max_speed_sfc'   }
    { 'DIAGS/hist_meas_speed_<CASE>.h5' '/s_max_speed'    '/s_max_speed'    '/s_max_speed_sfc'    }

    { 'DIAGS/hist_meas_speed_<CASE>.h5' '/ps_max_speed_t' '/ps_max_speed_t' '/ps_max_speed_t_sfc' }
    { 'DIAGS/hist_meas_speed_<CASE>.h5' '/s_max_speed_t'  '/s_max_speed_t'  '/s_max_speed_t_sfc'  }

    };
  Nsets = length(MeasList);

  Zstart = 0;     % layer to do the search for the max speed 
  Zend   = 1500;

  for icase = 1:Ncases
    Case = Cases{icase};

    fprintf('**********************************************************\n');
    fprintf('Generating max speed measurements for case: %s\n', Case);
    fprintf('\n');

    % Place all vars into one case specific output file.
    % If the file exists, remove it so that the HDF5 commands
    % can effectively re-create datasets.
    OutFile = sprintf('DIAGS/hist_meas_max_speed_%s.h5', Case);
    if (exist(OutFile, 'file') == 2)
      delete(OutFile);
    end

    for iset = 1:Nsets
      InFile      = regexprep(MeasList{iset}{1}, '<CASE>', Case);
      InVname     = MeasList{iset}{2};
      OutVname    = MeasList{iset}{3};
      OutSfcVname = MeasList{iset}{4};

      % Read in the variables
      fprintf('  Reading: %s (%s)\n', InFile, InVname);

      VAR = squeeze(h5read(InFile, InVname));
      X   = squeeze(h5read(InFile, '/x_coords'));
      Y   = squeeze(h5read(InFile, '/y_coords'));
      Z   = squeeze(h5read(InFile, '/z_coords'));
      T   = squeeze(h5read(InFile, '/t_coords'));

      Nx = length(X);
      Ny = length(Y);
      Nz = length(Z);
      Nt = length(T);

      Z1 = find(Z >= Zstart, 1, 'first');
      Z2 = find(Z <= Zend,   1, 'last');

      % If this is the first set, write out the coordinates into the output file
      % so that subsequent variables can have these coordinates attached.
      if (iset == 1)
        Xname = '/x_coords';
        Yname = '/y_coords';
        Zname = '/z_coords';
        Tname = '/t_coords';
  
        CreateDimensionsXyzt(OutFile, X, Y, Z, [ 1 ], Xname, Yname, Zname, Tname);
        % Add COARDS annotations
        NotateDimensionsXyzt(OutFile, Xname, Yname, Zname, Tname);
      end

      % Input fields are (r,z). Want to peel off surface layer (k == 2), and
      % the near surface layer that contains the speed max. Find this layer
      % by searching for the max value in the 2D (r,z) field. Do the max function
      % on the r dimension first so that the second max is looking along the z
      % dimension.
      [ MAX_VAL MAX_Z ] = max(max(VAR(:,Z1:Z2),[],1));

      VAR_SFC = squeeze(VAR(:,2));
      VAR_MAX = squeeze(VAR(:,MAX_Z));

      Vsize = Nx;
      DimOrder = { 'x' };

      % Collect output vars into one case specific file.
      fprintf('  Writing: %s (%s)\n', OutFile, OutVname);
      h5create(OutFile, OutVname, Vsize);
      h5write(OutFile, OutVname, VAR_MAX);

      fprintf('  Writing: %s (%s)\n', OutFile, OutSfcVname);
      h5create(OutFile, OutSfcVname, Vsize);
      h5write(OutFile, OutSfcVname, VAR_SFC);

      AttachDimensionsXyzt(OutFile, OutVname, DimOrder, Xname, Yname, Zname, Tname);
      AttachDimensionsXyzt(OutFile, OutSfcVname, DimOrder, Xname, Yname, Zname, Tname);

      fprintf('\n');
    end
  end
end
