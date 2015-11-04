function [ ] = GenTsavgHists()
% GenTsavgHists function to calculate temperature histograms

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
%    { 'HDF5/<CASE>/HDF5/lead_tempc25m_anom-<CASE>-AS-2006-08-20-120000-g3.h5' '/tempc' [ -5:.1:0 ] 'TsAveragedData/hist_lead_cpools_<CASE>.h5', '/cpools' }
%    { 'HDF5/<CASE>/HDF5/spath_tempc25m_anom-<CASE>-AS-2006-08-20-120000-g3.h5' '/tempc' [ -5:.1:0 ] 'TsAveragedData/hist_spath_cpools_<CASE>.h5', '/cpools' }
%    { 'HDF5/<CASE>/HDF5/smaxcp_tempc25m_anom-<CASE>-AS-2006-08-20-120000-g3.h5' '/tempc' [ -5:.1:0 ] 'TsAveragedData/hist_smaxcp_cpools_<CASE>.h5', '/cpools' }

    { 'HDF5/<CASE>/HDF5/core_tempc25m_anom-<CASE>-AS-2006-08-20-120000-g3.h5' '/tempc' [ -10:.1:0 ] 'TsAveragedData/hist_core_cpools_<CASE>.h5', '/cpools' }
    { 'HDF5/<CASE>/HDF5/rb_tempc25m_anom-<CASE>-AS-2006-08-20-120000-g3.h5' '/tempc' [ -10:.1:0 ] 'TsAveragedData/hist_rb_cpools_<CASE>.h5', '/cpools' }
    { 'HDF5/<CASE>/HDF5/env_tempc25m_anom-<CASE>-AS-2006-08-20-120000-g3.h5' '/tempc' [ -10:.1:0 ] 'TsAveragedData/hist_env_cpools_<CASE>.h5', '/cpools' }
    };
  Nfiles = length(FileList);

  for icase = 1:Ncases
    Case = CaseList{icase};

    fprintf('*****************************************************************\n');
    fprintf('Generating temperature histograms:\n');
    fprintf('  Case: %s\n', Case);
    fprintf('\n');

    for ifile = 1:Nfiles
      Ifile    = FileList{ifile}{1};
      Ivname   = FileList{ifile}{2};
      Bins     = FileList{ifile}{3};
      Ofile    = FileList{ifile}{4};
      Ovname   = FileList{ifile}{5};

      InFile  = regexprep(Ifile, '<CASE>', Case);
      OutFile = regexprep(Ofile, '<CASE>', Case);

      fprintf('  Reading: %s (%s)\n', InFile, Ivname);
      fprintf('\n');

      % HDATA will be (x,y,t)
      HDATA = squeeze(h5read(InFile, Ivname));
      X = squeeze(h5read(InFile, '/x_coords'));
      Y = squeeze(h5read(InFile, '/y_coords'));
      Z = squeeze(h5read(InFile, '/z_coords'));
      T = squeeze(h5read(InFile, '/t_coords'));

      Nx = length(X);
      Ny = length(Y);
      Nz = length(Z);
      Nt = length(T);

      Nb = length(Bins);

      % HDATA is (x,y,t) or (x,y,z,t)
      % nans will be ignored by histc
      % do histc across x dimension, then sum counts across y dimension
      HDATA = histc(HDATA, Bins, 1);
      HDATA = squeeze(sum(HDATA, 2));

      % Write out histogram data, mimic format of tsavg -> (x,y,z,t) where x is bins, y is dummy dimension
      fprintf('  Writing: %s (%s)\n', OutFile, Ovname)
      fprintf('\n');

      % get old file out of the way
      if (exist(OutFile, 'file') == 2)
        delete(OutFile);
      end

      % HDATA is (b,t) or (b,z,t) now, need to expand into (x,y,z,t);
      HDATA = reshape(HDATA, [ Nb 1 Nz Nt ]);
      DimOrder = { 'x' 'y' 'z' 't' };
      Vsize = [ Nb 1 Nz Nt ];

      h5create(OutFile, Ovname, Vsize);
      h5write (OutFile, Ovname, HDATA);

      Xname = '/x_coords';
      Yname = '/y_coords';
      Zname = '/z_coords';
      Tname = '/t_coords';

      CreateDimensionsXyzt(OutFile, Bins, [1], Z, T, Xname, Yname, Zname, Tname);
      % Add COARDS annotations
      NotateDimensionsXyzt(OutFile, Xname, Yname, Zname, Tname);

      % Attach dimensions
      AttachDimensionsXyzt(OutFile, Ovname, DimOrder, Xname, Yname, Zname, Tname);
    end
  end
end 
