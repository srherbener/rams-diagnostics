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
    { 'HDF5/<CASE>/HDF5/u-<CASE>-AS-2006-08-20-120000-g3.h5' 'HDF5/<CASE>/HDF5/u_lite-<CASE>-AS-2006-08-20-120000-g3.h5' '/u' 10 10 1 1 }
    { 'HDF5/<CASE>/HDF5/v-<CASE>-AS-2006-08-20-120000-g3.h5' 'HDF5/<CASE>/HDF5/v_lite-<CASE>-AS-2006-08-20-120000-g3.h5' '/v' 10 10 1 1 }
    };

  Nfiles = length(FileList);

  for icase = 1:Ncases
    Case = CaseList{icase};

    fprintf('*****************************************************************\n');
    fprintf('Reducing data for case: %s\n', Case);
    fprintf('\n');


    for ifile = 1:Nfiles
      InFile  = regexprep(FileList{ifile}{1}, '<CASE>', Case);
      OutFile = regexprep(FileList{ifile}{2}, '<CASE>', Case);
      InVname = FileList{ifile}{3};
      Xinc    = FileList{ifile}{4};
      Yinc    = FileList{ifile}{5};
      Zinc    = FileList{ifile}{6};
      Tinc    = FileList{ifile}{7};

      ReduceData(InFile, OutFile, InVname, Xinc, Yinc, Zinc, Tinc);

    end % files
  end % cases
end % function
