function [ ] = GenAzHistSumsCounts()

  % make sure output directory exists
  Ddir = 'DIAGS';  % coordinate this with output file names below
  if (exist(Ddir, 'dir') ~= 7)
    mkdir(Ddir);
  end

  % list of simulation cases
  CaseList = {
    'TSD_SAL_DUST'
%    'TSD_SAL_NODUST'
%    'TSD_NONSAL_DUST'
%    'TSD_NONSAL_NODUST'
    };
  Ncases = length(CaseList);


  % measurement sets
  MeasSets = {

    % Dust
    {
      'Dust Azavg'
      {

        % Region of storm (500 km radius)
        { 'AzAveragedData/hist_all_d1_mass_<CASE>.h5'     '/d1_mass'     '/all_sum_d1_mass'     } 
        { 'AzAveragedData/hist_all_d2_mass_<CASE>.h5'     '/d2_mass'     '/all_sum_d2_mass'     } 
        { 'AzAveragedData/hist_all_dust_mass_<CASE>.h5'   '/dust_mass'   '/all_sum_dust_mass'   } 
        { 'AzAveragedData/hist_all_dust_hydro_<CASE>.h5'  '/dust_hydro'  '/all_sum_dust_hydro'  }
        { 'AzAveragedData/hist_all_tracer_mass_<CASE>.h5' '/tracer_mass' '/all_sum_tracer_mass' }
        { 'AzAveragedData/hist_all_dust_sfc_<CASE>.h5'    '/dust_sfc'    '/all_sum_dust_sfc'    }
      }
      'DIAGS/hist_sums_az_dust_<CASE>.h5'
    }

    % CCN
    {
      'CCN Azavg'
      {
        % Region of storm (500 km radius)
        { 'AzAveragedData/hist_all_ccn_mass_<CASE>.h5'    '/ccn_mass'    '/all_sum_ccn_mass' }
      }
      'DIAGS/hist_sums_az_ccn_<CASE>.h5'
    }

    % Regen
    {
      'Regen Azavg'
      {
        % Region of storm (500 km radius)
        { 'AzAveragedData/hist_all_ra_mass_<CASE>.h5'     '/ra_mass'     '/all_sum_ra_mass' }
      }
      'DIAGS/hist_sums_az_ra_<CASE>.h5'
    }

    % All aerosols
    {
      'Aero Azavg'
      {
        % Region of storm (500 km radius)
        { 'AzAveragedData/hist_all_aero_mass_<CASE>.h5'   '/aero_mass'   '/all_sum_aero_mass' }
      }
      'DIAGS/hist_sums_az_aero_<CASE>.h5'
    }

    };

  Nsets = length(MeasSets);

  for icase = 1:Ncases
    Case = CaseList{icase};

    fprintf('*****************************************************************\n');
    fprintf('Generating histogram sums and counts:\n');
    fprintf('  Case: %s\n', Case);
    fprintf('\n');


    for iset = 1:Nsets
      MeasName  = MeasSets{iset}{1};
      MeasList  = MeasSets{iset}{2};
      OutFname  = MeasSets{iset}{3};

      fprintf('    Measurement set: %s\n', MeasName);
      fprintf('\n');

      % Put all measurements into one file per case
      % If the file exists, remove it so that the HDF5 commands
      % can effectively re-create datasets.
      OutFile = regexprep(OutFname, '<CASE>', Case);
      if (exist(OutFile, 'file') == 2)
        delete(OutFile);
      end
  
      icount = 0;
      Nmeas = length(MeasList);
      for imeas = 1:Nmeas
        InFname = MeasList{imeas}{1};
        InVname     = MeasList{imeas}{2};
        OutVname  = MeasList{imeas}{3};

        % skip this profile set if doing dust and on a NODUST case
        if ((~isempty(regexp(Case, 'NODUST'))) && ...
            ((~isempty(regexp(InVname, '/d[12]_num'))) || ...
             (~isempty(regexp(InVname, '/d[12]_mass'))) || ...
             (~isempty(regexp(InVname, '/dust_'))) || ...
             (~isempty(regexp(InVname, '/dustifn_'))) || ...
             (~isempty(regexp(InVname, '/tracer[12]'))) || ...
             (~isempty(regexp(InVname, '/trdust[12]_diff')))))
          continue
        else
          icount = icount + 1;
        end

        InFile = regexprep(InFname, '<CASE>', Case);
        fprintf('      Reading: %s (%s)\n', InFile, InVname);

        HDATA = squeeze(h5read(InFile, InVname));
        X     = squeeze(h5read(InFile, '/x_coords'));
        Y     = squeeze(h5read(InFile, '/y_coords'));
        Z     = squeeze(h5read(InFile, '/z_coords'));
        T     = squeeze(h5read(InFile, '/t_coords'));

        Nx = length(X);
        Ny = length(Y);
        Nz = length(Z);
        Nt = length(T);

        % HDATA is of the form:
        %    2D field: (x,y,t)
        %    3D field: (x,y,z,t)
        %
        % where:
        %   x -> radial bands
        %   y -> histogram bins
        %   z -> height
        %   t -> time
        FieldIs3D = ndims(HDATA) == 4;

        % Reduce x dimension (radial bands) by summing up the histogram counts
        HDATA = squeeze(sum(HDATA,1));

        % HDATA is now either (y,z,t) or (y,t)
        %
        % Reduce the y dimension (histogram counts) by
        %   SUMS = sum of each count * bin value
        %   COUNTS = sum of each count
        %
        % Bins (Y) are the edges of the bins, so use the average values between the edges
        % for "bin value". The last entry in Y was used as an exact match during the histogram
        % counting so use it for the final entry in BIN_VALS.
        BIN_VALS = (Y(1:end-1) + Y(2:end)) .* 0.5;
        BIN_VALS(Ny) = Y(Ny);

        if (FieldIs3D)
          BIN_VALS = repmat(BIN_VALS, [ 1 Nz Nt ]);
        else
          BIN_VALS = repmat(BIN_VALS, [ 1 Nt ]);
        end

        % Create two output variables:
        %   SUMS: contains sums across horizontal domain
        %   COUNTS: contains total of each bin count
        SUMS = squeeze(sum((HDATA .* BIN_VALS), 1));
        COUNTS = squeeze(sum(HDATA, 1));

        % If first measurement, then write out coordinates for later
        % use in attaching vars to them. Write out the original coordinates.
        % The original coordinates with their original lengths will always be the
        % appropriate vectors to use for attaching to variables.
        if (icount == 1)
          Xname = '/x_coords';
          Yname = '/y_coords';
          Zname = '/z_coords';
          Tname = '/t_coords';
  
          CreateDimensionsXyzt(OutFile, X, Y, Z, T, Xname, Yname, Zname, Tname);
          % Add COARDS annotations
          NotateDimensionsXyzt(OutFile, Xname, Yname, Zname, Tname);
        end

        % SUMS and COUNTS:
        %   2D field: (t)
        %   3d field: (z,t)
        if (FieldIs3D)
          OutSize = [ Nz Nt ];
          DimOrder = { 'z' 't' };
        else
          OutSize = Nt;
          DimOrder = { 't' };
        end

        OutVnameCounts = sprintf('%s_counts', OutVname);
        fprintf('      Writing: %s (%s)\n', OutFile, OutVname)
        fprintf('      Writing: %s (%s)\n', OutFile, OutVnameCounts)

        h5create(OutFile, OutVname, OutSize);
        h5write(OutFile, OutVname, SUMS);

        h5create(OutFile, OutVnameCounts, OutSize);
        h5write(OutFile, OutVnameCounts, COUNTS);

        % Attach dimensions
        AttachDimensionsXyzt(OutFile, OutVname, DimOrder, Xname, Yname, Zname, Tname);
        AttachDimensionsXyzt(OutFile, OutVnameCounts, DimOrder, Xname, Yname, Zname, Tname);

        fprintf('\n');
      end % measurements
    end % sets
  end % cases
end % function
