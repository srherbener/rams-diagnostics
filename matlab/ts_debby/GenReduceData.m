function [ ] = GenReduceData()

  % list of simulation cases
  CaseList = {
    'TSD_SAL_DUST'
    'TSD_SAL_NODUST'
    'TSD_NONSAL_DUST'
    'TSD_NONSAL_NODUST'
    };
  Ncases = length(CaseList);

  FileList = {
    % Horizontal winds
    { 'HDF5/<CASE>/HDF5/u-<CASE>-AS-2006-08-20-120000-g3.h5' 'HDF5/<CASE>/HDF5/u_lite-<CASE>-AS-2006-08-20-120000-g3.h5' '/u' 10 10 1 1 'm/s' 'u' 't z y x' }
    { 'HDF5/<CASE>/HDF5/v-<CASE>-AS-2006-08-20-120000-g3.h5' 'HDF5/<CASE>/HDF5/v_lite-<CASE>-AS-2006-08-20-120000-g3.h5' '/v' 10 10 1 1 'm/s' 'u' 't z y x' }

    % sea_pressure
    { 'HDF5/<CASE>/HDF5/sea_press-<CASE>-AS-2006-08-20-120000-g3.h5' 'HDF5/<CASE>/HDF5/sea_press_lite-<CASE>-AS-2006-08-20-120000-g3.h5' '/sea_press' 10 10 1 1 'mb' 'sea-level-pressure' 't y x' }

    % topography
    { 'HDF5/<CASE>/HDF5/topo-<CASE>-AS-2006-08-20-120000-g3.h5' 'HDF5/<CASE>/HDF5/topo_lite-<CASE>-AS-2006-08-20-120000-g3.h5' '/topo' 10 10 1 1 'm' 'topo' 't y x' }
    };

  Nfiles = length(FileList);

  for icase = 1:Ncases
    Case = CaseList{icase};

    fprintf('*****************************************************************\n');
    fprintf('Reducing data for case: %s\n', Case);
    fprintf('\n');


    for ifile = 1:Nfiles
      InFile   = regexprep(FileList{ifile}{1}, '<CASE>', Case);
      OutFile  = regexprep(FileList{ifile}{2}, '<CASE>', Case);
      InVname  = FileList{ifile}{3};
      Xinc     = FileList{ifile}{4};
      Yinc     = FileList{ifile}{5};
      Zinc     = FileList{ifile}{6};
      Tinc     = FileList{ifile}{7};
      Units    = FileList{ifile}{8};
      LongName = FileList{ifile}{9};
      DimNames = FileList{ifile}{10};

      ReduceData(InFile, OutFile, InVname, Xinc, Yinc, Zinc, Tinc, Units, LongName, DimNames);

    end % files
  end % cases
end % function
