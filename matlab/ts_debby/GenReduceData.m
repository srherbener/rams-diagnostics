function [ ] = GenReduceData()

  % list of simulation cases
  CaseList = {
%    'TSD_SAL_DUST'
%    'TSD_SAL_NODUST'
    'TSD_NONSAL_DUST'
%    'TSD_NONSAL_NODUST'
    };
  Ncases = length(CaseList);

  FileList = {
%    % Dust
%    { 'HDF5/<CASE>/HDF5/dust_mass-<CASE>-AS-2006-08-20-120000-g3.h5' 'HDF5/<CASE>/HDF5/dust_mass_lite-<CASE>-AS-2006-08-20-120000-g3.h5' '/dust_mass' 10 10 1 1 'C' 'dust mass' 't z y x' }
%
%    % Cloud
%    { 'HDF5/<CASE>/HDF5/cloud-<CASE>-AS-2006-08-20-120000-g3.h5' 'HDF5/<CASE>/HDF5/cloud_lite-<CASE>-AS-2006-08-20-120000-g3.h5' '/cloud' 10 10 1 1 'C' 'cloud' 't z y x' }
%
%    % Total Condensate
%    { 'HDF5/<CASE>/HDF5/total_cond-<CASE>-AS-2006-08-20-120000-g3.h5' 'HDF5/<CASE>/HDF5/total_cond_lite-<CASE>-AS-2006-08-20-120000-g3.h5' '/total_cond' 10 10 1 1 'C' 'total condendsate' 't z y x' }
%
%%    % Temp
%    { 'HDF5/<CASE>/HDF5/tempc-<CASE>-AS-2006-08-20-120000-g3.h5' 'HDF5/<CASE>/HDF5/tempc_lite-<CASE>-AS-2006-08-20-120000-g3.h5' '/tempc' 10 10 1 1 'C' 'tempc' 't z y x' }
%    { 'HDF5/<CASE>/HDF5/theta-<CASE>-AS-2006-08-20-120000-g3.h5' 'HDF5/<CASE>/HDF5/theta_lite-<CASE>-AS-2006-08-20-120000-g3.h5' '/theta' 10 10 1 1 'C' 'theta' 't z y x' }
%    { 'HDF5/<CASE>/HDF5/theta_v-<CASE>-AS-2006-08-20-120000-g3.h5' 'HDF5/<CASE>/HDF5/theta_v_lite-<CASE>-AS-2006-08-20-120000-g3.h5' '/theta_v' 10 10 1 1 'C' 'theta_v' 't z y x' }
%
%    % Vapor
%    { 'HDF5/<CASE>/HDF5/vapor-<CASE>-AS-2006-08-20-120000-g3.h5' 'HDF5/<CASE>/HDF5/vapor_lite-<CASE>-AS-2006-08-20-120000-g3.h5' '/vapor' 10 10 1 1 'C' 'vapor' 't z y x' }
%
%    % Horizontal winds
%    { 'HDF5/<CASE>/HDF5/u-<CASE>-AS-2006-08-20-120000-g3.h5' 'HDF5/<CASE>/HDF5/u_lite-<CASE>-AS-2006-08-20-120000-g3.h5' '/u' 10 10 1 1 'm/s' 'u' 't z y x' }
%    { 'HDF5/<CASE>/HDF5/v-<CASE>-AS-2006-08-20-120000-g3.h5' 'HDF5/<CASE>/HDF5/v_lite-<CASE>-AS-2006-08-20-120000-g3.h5' '/v' 10 10 1 1 'm/s' 'u' 't z y x' }
%
%    % sea_pressure
    { 'HDF5/<CASE>/HDF5/sea_press-<CASE>-AS-2006-08-20-120000-g3.h5' 'HDF5/<CASE>/HDF5/sea_press_lite-<CASE>-AS-2006-08-20-120000-g3.h5' '/sea_press' 10 10 1 1 'mb' 'sea-level-pressure' 't y x' }

%    % topography
    { 'HDF5/<CASE>/HDF5/topo-<CASE>-AS-2006-08-20-120000-g3.h5' 'HDF5/<CASE>/HDF5/topo_lite-<CASE>-AS-2006-08-20-120000-g3.h5' '/topo' 10 10 1 1 'm' 'topo' 't y x' }

%    % Horizontal winds - pressure surfaces
    { 'HDF5/<CASE>/HDF5/u-<CASE>-AP-2006-08-20-120000-g3.h5' 'HDF5/<CASE>/HDF5/u_lite-<CASE>-AP-2006-08-20-120000-g3.h5' '/u' 10 10 1 1 'm/s' 'u' 't z y x' }
    { 'HDF5/<CASE>/HDF5/v-<CASE>-AP-2006-08-20-120000-g3.h5' 'HDF5/<CASE>/HDF5/v_lite-<CASE>-AP-2006-08-20-120000-g3.h5' '/v' 10 10 1 1 'm/s' 'u' 't z y x' }

%    % Temp - pressure surfaces
    { 'HDF5/<CASE>/HDF5/tempc-<CASE>-AP-2006-08-20-120000-g3.h5' 'HDF5/<CASE>/HDF5/tempc_lite-<CASE>-AP-2006-08-20-120000-g3.h5' '/tempc' 10 10 1 1 'C' 'tempc' 't z y x' }

    % Vapor - pressure surfaces
    { 'HDF5/<CASE>/HDF5/vapor-<CASE>-AP-2006-08-20-120000-g3.h5' 'HDF5/<CASE>/HDF5/vapor_lite-<CASE>-AP-2006-08-20-120000-g3.h5' '/vapor' 10 10 1 1 'C' 'vapor' 't z y x' }

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

      switch InVname 
        case { '/dust_mass' }
          if (regexp(Case, '_NODUST'))
            fprintf('  WARNING: attempting to process dust variable in non-dust case: skipping this case\n');
            continue;
          end
      end

      ReduceData(InFile, OutFile, InVname, Xinc, Yinc, Zinc, Tinc, Units, LongName, DimNames);

    end % files
  end % cases
end % function
