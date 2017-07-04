function [ ] = GenSalLocs()
% This function will generate SAL location files for a given TS Debby simulation

  LonRange = [ -50  5 ]; % only look in box defined by LonRange, LatRange, and Zrange
  LatRange = [  10 32 ];
  Zrange   = [   1  5 ]; 
  
  ColCount = 5; % Have to have at least 5 "hits" per column according to the test
                % for less than RhLimit to count the column as part of the SAL.
  
  FileList = {
    % orig sims
    { 'HDF5/TSD_SAL_DUST/HDF5/relhum-TSD_SAL_DUST-AS-2006-08-20-120000-g1.h5' '/relhum' 1 45 'HDF5/namma_profiles/sal_loc-TSD_DRY_DUST-AS-2006-08-20-120000-g1.h5' }
    { 'HDF5/TSD_SAL_DUST/HDF5/relhum-TSD_SAL_DUST-AS-2006-08-20-120000-g2.h5' '/relhum' 1 45 'HDF5/namma_profiles/sal_loc-TSD_DRY_DUST-AS-2006-08-20-120000-g2.h5' }
    { 'HDF5/TSD_SAL_DUST/HDF5/relhum-TSD_SAL_DUST-AS-2006-08-20-120000-g3.h5' '/relhum' 1 45 'HDF5/namma_profiles/sal_loc-TSD_DRY_DUST-AS-2006-08-20-120000-g3.h5' }
  
    % new 3 grid sim
    { 'HDF5/MID_LEVEL/TSD_SPIN_UP/HDF5/relhum-TSD_SPIN_UP-AS-2006-08-20-120000-g1.h5' '/relhum' 1 48 'HDF5/namma_profiles/sal_loc-TSD_SPIN_UP-AS-2006-08-20-120000-g1.h5' }
    { 'HDF5/MID_LEVEL/TSD_SPIN_UP/HDF5/relhum-TSD_SPIN_UP-AS-2006-08-20-120000-g2.h5' '/relhum' 1 48 'HDF5/namma_profiles/sal_loc-TSD_SPIN_UP-AS-2006-08-20-120000-g2.h5' }
    { 'HDF5/MID_LEVEL/TSD_SPIN_UP/HDF5/relhum-TSD_SPIN_UP-AS-2006-08-20-120000-g3.h5' '/relhum' 1 48 'HDF5/namma_profiles/sal_loc-TSD_SPIN_UP-AS-2006-08-20-120000-g3.h5' }
  
    % new 1 grid sim
    { 'HDF5/MID_LEVEL/TSD_SPIN_UP_1G/HDF5/relhum-TSD_SPIN_UP_1G-AS-2006-08-21-180000-g1.h5' '/relhum' 13 47 'HDF5/namma_profiles/sal_loc-TSD_SPIN_UP_1G-AS-2006-08-21-180000-g1.h5' }
    };
  Nfiles = length(FileList);
  
  fprintf('Generating SAL location files:\n');
  for ifile = 1:Nfiles
    InFile   = FileList{ifile}{1};
    InVname  = FileList{ifile}{2};
    TimeStep = FileList{ifile}{3};
    RhLimit  = FileList{ifile}{4};
    OutFile  = FileList{ifile}{5};
    
    LocateSal(InFile, InVname, TimeStep, RhLimit, LonRange, LatRange, Zrange, ColCount, OutFile);

    fprintf('\n');
  end
end
