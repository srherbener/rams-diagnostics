function [ ] = GenAccumSfcDust()
% GenAccumSfcDust function to calculate total accum dust from a histogram extraction

  % Cases
  CaseList = {
    'TSD_SAL_DUST'
    'TSD_NONSAL_DUST'
    };
  Ncases = length(CaseList);

  % Description of measurements
  FileList = {

    { 'AzAveragedData/hist_all_dust_sfc_<CASE>.h5'   '/dust_sfc' 'DIAGS/hist_meas_az_dust_sfc_<CASE>.h5' '/all_dust_sfc'   'ug/m3' 'surface_accum_dust' 't' }
    { 'TsAveragedData/hist_spath_dust_sfc_<CASE>.h5' '/dust_sfc' 'DIAGS/hist_meas_ts_dust_sfc_<CASE>.h5' '/spath_dust_sfc' 'ug/m3' 'surface_accum_dust' 't' }

    };
  Nfiles = length(FileList);

  for icase = 1:Ncases
    Case = CaseList{icase};

    fprintf('*****************************************************************\n');
    fprintf('Generating surface accumulated dust amounts:\n');
    fprintf('  Case: %s\n', Case);
    fprintf('\n');

    for ifile = 1:Nfiles
      Dfile    = FileList{ifile}{1};
      Dvname   = FileList{ifile}{2};
      Ofile    = FileList{ifile}{3};
      Ovname   = FileList{ifile}{4};
      Units    = FileList{ifile}{5};
      LongName = FileList{ifile}{6};
      DimNames = FileList{ifile}{7};

      InFile  = regexprep(Dfile, '<CASE>', Case);
      OutFile = regexprep(Ofile, '<CASE>', Case);

      fprintf('  Reading: %s (%s)\n', InFile, Dvname);
      fprintf('\n');

      % DDATA for AzAveragedData is: (x,y,t)
      %   x - radius
      %   y - bin edges
      %   t - time
      %
      % DDATA for TsAveragedData is: (y,t)
      %   y - bin edges
      %   t - time
      DDATA = squeeze(h5read(InFile, Dvname));
      X = squeeze(h5read(InFile, '/x_coords'));
      Y = squeeze(h5read(InFile, '/y_coords'));
      Z = squeeze(h5read(InFile, '/z_coords'));
      T = squeeze(h5read(InFile, '/t_coords'));

      Nx = length(X);
      Ny = length(Y);
      Nz = length(Z);
      Nt = length(T);

      % Bin values are the edges of the bins. Take the average between successive edges to
      % represent the value for each bin. On the last bin, just use the value in X(Nx).
      if (~isempty(regexp(InFile, 'AzAveragedData')))
        BinVals = (Y(1:end-1) + Y(2:end)) .* 0.5;
        BinVals(Ny) = Y(Ny);
      else
        BinVals = (X(1:end-1) + X(2:end)) .* 0.5;
        BinVals(Nx) = X(Nx);
      end
      BVALS = repmat(BinVals, [ 1 Nt ]);

      % Sum up the histogram counts across the x (radius) domain, then convert the resultant histogram
      % to an accumulated value.
      % Only do the first sum (reduction of x dimension) if AzAveragedData
      if (~isempty(regexp(InFile, 'AzAveragedData')))
        DDATA = squeeze(sum(DDATA,1));
      end
      ADUST = squeeze(sum((DDATA .* BVALS),1));

      % Write out anomaly data
      fprintf('  Writing: %s (%s)\n', OutFile, Ovname)
      fprintf('\n');

      % get old file out of the way
      if (exist(OutFile, 'file') == 2)
        delete(OutFile);
      end

      % Write out each measurement (time series)
      % Need to write out 4D var (add in the trivial Z dimension)
      DimOrder = { 't' };
      Vsize = Nt;

      h5create(OutFile, Ovname, Vsize, 'DataType', 'single');
      h5write (OutFile, Ovname, ADUST);
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
