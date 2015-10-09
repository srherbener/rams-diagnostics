function [ ] = GenTempAnomalies()
% GenTempAnomalies function to calculate temperature anomalies

  % Cases
  CaseList = {
    'TSD_SAL_DUST'
    'TSD_SAL_NODUST'
    'TSD_NONSAL_DUST'
    'TSD_NONSAL_NODUST'
    };
  Ncases = length(CaseList);

  % Description of measurements
  FileList = {
    { 'HDF5/<CASE>/HDF5/lead_tempc25m-<CASE>-AS-2006-08-20-120000-g3.h5' '/tempc' 'HDF5/<CASE>/HDF5/lead_tempc25m_anom-<CASE>-AS-2006-08-20-120000-g3.h5' '/tempc' 'C' 'temperature' 't z y x' }
    };
  Nfiles = length(FileList);

  for icase = 1:Ncases
    Case = CaseList{icase};

    fprintf('*****************************************************************\n');
    fprintf('Generating temperature anomalies:\n');
    fprintf('  Case: %s\n', Case);
    fprintf('\n');

    for ifile = 1:Nfiles
      Tfile    = FileList{ifile}{1};
      Tvname   = FileList{ifile}{2};
      Ofile    = FileList{ifile}{3};
      Ovname   = FileList{ifile}{4};
      Units    = FileList{ifile}{5};
      LongName = FileList{ifile}{6};
      DimNames = FileList{ifile}{7};

      InFile  = regexprep(Tfile, '<CASE>', Case);
      OutFile = regexprep(Ofile, '<CASE>', Case);

      fprintf('  Reading: %s (%s)\n', InFile, Tvname);
      fprintf('\n');

      % TDATA will be 2D field, ie. (x,y,t)
      TDATA = squeeze(h5read(InFile, Tvname));
      X = squeeze(h5read(InFile, '/x_coords'));
      Y = squeeze(h5read(InFile, '/y_coords'));
      Z = squeeze(h5read(InFile, '/z_coords'));
      T = squeeze(h5read(InFile, '/t_coords'));

      Nx = length(X);
      Ny = length(Y);
      Nz = length(Z);
      Nt = length(T);

      % Calculate mean, leaving out zero values, then subract mean from non-zero input values
      % TDATA is (x,y,t);
      TEMP = TDATA;
      TEMP(TEMP == 0) = nan;
      AVG = nanmean(nanmean(TEMP,2),1); % time series of mean, [ 1 1 Nt ]
      AVG = repmat(AVG, [ Nx Ny 1 ]); % repeat time series of means over horiz domain

      % TEMP has nans in the regions that we don't want the calculation to take place
      % so allow these to remain in output.
      TEMP = TEMP - AVG;

      % Write out anomaly data
      fprintf('  Writing: %s (%s)\n', OutFile, Ovname)
      fprintf('\n');

      % get old file out of the way
      if (exist(OutFile, 'file') == 2)
        delete(OutFile);
      end

      % Write out each measurement (time series)
      % Need to write out 4D var (add in the trivial Z dimension)
      TEMP = reshape(TEMP, [ Nx Ny 1 Nt ]);
      DimOrder = { 'x' 'y' 'z' 't' };
      Vsize = [ Nx Ny 1 Nt ];

      h5create(OutFile, Ovname, Vsize, 'ChunkSize', [ Nx Ny 1 1 ], 'Deflate', 1, 'Shuffle', 1);
      h5write (OutFile, Ovname, TEMP);
      NotateVariableXyzt(OutFile, Ovname, Units, LongName, DimNames);

      Xname = '/x_coords';
      Yname = '/y_coords';
      Zname = '/z_coords';
      Tname = '/t_coords';

      CreateDimensionsXyzt(OutFile, X, Y, Z, T, Xname, Yname, Zname, Tname);
      % Add COARDS annotations
      NotateDimensionsXyzt(OutFile, Xname, Yname, Zname, Tname);

      % Attach dimensions
      AttachDimensionsXyzt(OutFile, Ovname, DimOrder, Xname, Yname, Zname, Tname);
    end
  end
end 
