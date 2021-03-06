function [ ] = GenStormDustFields()
% GenStormDustFields function to select dust fields only within a 500 km radius of the storm center.

  % Cases
  CaseList = {
    'TSD_SAL_DUST'
    'TSD_NONSAL_DUST'
    };
  Ncases = length(CaseList);

  % Description of dust fields
  FileList = {
    { 'FILTERS/all_500_<CASE>.h5' '/filter' 'HDF5/<CASE>/HDF5/d1_mass-<CASE>-AS-2006-08-20-120000-g3.h5' '/dust1_mass' 'VAPOR/HDF5/storm_d1_mass_<CASE>.h5', '/d1_mass' }
    { 'FILTERS/all_500_<CASE>.h5' '/filter' 'HDF5/<CASE>/HDF5/d2_mass-<CASE>-AS-2006-08-20-120000-g3.h5' '/dust2_mass' 'VAPOR/HDF5/storm_d2_mass_<CASE>.h5', '/d2_mass' }

    { 'FILTERS/all_500_<CASE>.h5' '/filter' 'HDF5/<CASE>/HDF5/d1_num-<CASE>-AS-2006-08-20-120000-g3.h5' '/dust1_concen' 'VAPOR/HDF5/storm_d1_num_<CASE>.h5', '/d1_num' }
    { 'FILTERS/all_500_<CASE>.h5' '/filter' 'HDF5/<CASE>/HDF5/d2_num-<CASE>-AS-2006-08-20-120000-g3.h5' '/dust2_concen' 'VAPOR/HDF5/storm_d2_num_<CASE>.h5', '/d2_num' }

    { 'FILTERS/all_500_<CASE>.h5' '/filter' 'HDF5/<CASE>/HDF5/trdust1_diff-<CASE>-AS-2006-08-20-120000-g3.h5' '/trdust1_diff' 'VAPOR/HDF5/storm_trdust1_diff_<CASE>.h5', '/trdust1_diff' }
    { 'FILTERS/all_500_<CASE>.h5' '/filter' 'HDF5/<CASE>/HDF5/trdust2_diff-<CASE>-AS-2006-08-20-120000-g3.h5' '/trdust2_diff' 'VAPOR/HDF5/storm_trdust2_diff_<CASE>.h5', '/trdust2_diff' }

    { 'FILTERS/all_500_<CASE>.h5' '/filter' 'HDF5/<CASE>/HDF5/dust_cloud-<CASE>-AS-2006-08-20-120000-g3.h5' '/dust_cloud_mass' 'VAPOR/HDF5/storm_dust_cloud_<CASE>.h5', '/dust_cloud' }
    { 'FILTERS/all_500_<CASE>.h5' '/filter' 'HDF5/<CASE>/HDF5/dust_rain-<CASE>-AS-2006-08-20-120000-g3.h5'  '/dust_rain_mass'  'VAPOR/HDF5/storm_dust_rain_<CASE>.h5',  '/dust_rain'  }

    { 'FILTERS/all_500_<CASE>.h5' '/filter' 'HDF5/<CASE>/HDF5/dust_pris-<CASE>-AS-2006-08-20-120000-g3.h5'  '/dust_pris_mass'  'VAPOR/HDF5/storm_dust_pris_<CASE>.h5',  '/dust_pris'  }
    { 'FILTERS/all_500_<CASE>.h5' '/filter' 'HDF5/<CASE>/HDF5/dust_snow-<CASE>-AS-2006-08-20-120000-g3.h5'  '/dust_snow_mass'  'VAPOR/HDF5/storm_dust_snow_<CASE>.h5',  '/dust_snow'  }
    { 'FILTERS/all_500_<CASE>.h5' '/filter' 'HDF5/<CASE>/HDF5/dust_aggr-<CASE>-AS-2006-08-20-120000-g3.h5'  '/dust_aggr_mass'  'VAPOR/HDF5/storm_dust_aggr_<CASE>.h5',  '/dust_aggr'  }
    { 'FILTERS/all_500_<CASE>.h5' '/filter' 'HDF5/<CASE>/HDF5/dust_graup-<CASE>-AS-2006-08-20-120000-g3.h5'  '/dust_grau_mass'  'VAPOR/HDF5/storm_dust_graup_<CASE>.h5',  '/dust_graup'  }
    { 'FILTERS/all_500_<CASE>.h5' '/filter' 'HDF5/<CASE>/HDF5/dust_hail-<CASE>-AS-2006-08-20-120000-g3.h5'  '/dust_hail_mass'  'VAPOR/HDF5/storm_dust_hail_<CASE>.h5',  '/dust_hail'  }
  };
  Nfiles = length(FileList);

  Xname = '/x_coords';
  Yname = '/y_coords';
  Zname = '/z_coords';
  Tname = '/t_coords';

  for icase = 1:Ncases
    Case = CaseList{icase};

    fprintf('*****************************************************************\n');
    fprintf('Generating storm dust fields:\n');
    fprintf('  Case: %s\n', Case);
    fprintf('\n');

    for ifile = 1:Nfiles
      Ffile    = FileList{ifile}{1};
      Fvname   = FileList{ifile}{2};
      Dfile    = FileList{ifile}{3};
      Dvname   = FileList{ifile}{4};
      Ofile    = FileList{ifile}{5};
      Ovname   = FileList{ifile}{6};

      FilterFile = regexprep(Ffile, '<CASE>', Case);
      DustFile   = regexprep(Dfile, '<CASE>', Case);
      OutFile    = regexprep(Ofile, '<CASE>', Case);

      fprintf('  Reading: %s (%s)\n', FilterFile, Fvname);
      fprintf('  Reading: %s (%s)\n', DustFile, Dvname);
      fprintf('  Writing: %s (%s)\n', OutFile, Ovname);
      fprintf('\n');

      X = squeeze(h5read(DustFile, Xname));
      Y = squeeze(h5read(DustFile, Yname));
      Z = squeeze(h5read(DustFile, Zname));
      T = squeeze(h5read(DustFile, Tname));

      Nx = length(X);
      Ny = length(Y);
      Nz = length(Z);
      Nt = length(T);

      % Get old output file out of the way, and write coordinate
      % values into the new file.
      if (exist(OutFile, 'file') == 2)
        delete(OutFile);
      end

      CreateDimensionsXyzt(OutFile, X, Y, Z, T, Xname, Yname, Zname, Tname);
      % Add COARDS annotations
      NotateDimensionsXyzt(OutFile, Xname, Yname, Zname, Tname);

      % Create output dataset, make the t dimension expandable 
      DimOrder = { 'x' 'y' 'z' 't' };
      Vsize = [ Nx Ny Nz Inf ];
      Csize = [ Nx Ny 1 1 ];
      h5create(OutFile, Ovname, Vsize, 'DataType', 'single', 'ChunkSize', Csize, 'Deflate', 6, 'Shuffle', true);

      % Do one time step at a time to help keep memory requirements lower.
      %   Filter is (x,y,t)
      %   Data is (x,y,z,t)
      for it = 1:Nt
        Start = [ 1 1 1 it ];
        Fcount = [ Nx Ny  1 1 ];
        Dcount = [ Nx Ny Nz 1 ];

        FILTER = squeeze(h5read(FilterFile, Fvname, Start, Fcount));
        DUST   = squeeze(h5read(DustFile,   Dvname, Start, Dcount));

        % Replicate levels in the filter so it matches up with the data
        FILTER = repmat(FILTER, [ 1 1 Nz ]);

        % Select dust only where filter has 1's (filter entries are either 1's or 0's)
        DUST = DUST .* FILTER;

        % Write DUST to output file. 
        h5write (OutFile, Ovname, DUST, Start, Dcount);

        if (mod(it, 10) == 0)
          fprintf('  Completed timestep: %d\n', it);
        end
      end
      
      % Attach dimensions
      AttachDimensionsXyzt(OutFile, Ovname, DimOrder, Xname, Yname, Zname, Tname);
    end
  end
end 
