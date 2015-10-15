function [ ] = GenStormDustFields()
% GenStormDustFields function to select dust fields only within a 500 km radius of the storm center.

  % Cases
  CaseList = {
%    'TSD_SAL_DUST'
    'TSD_NONSAL_DUST'
    };
  Ncases = length(CaseList);

  % Description of dust fields
  FileList = {
%    { 'FILTERS/all_500_TSD_SAL_DUST.h5' '/filter' 'HDF5/<CASE>/HDF5/d1_mass-<CASE>-AS-2006-08-20-120000-g3.h5' '/dust1_mass' 'VAPOR/HDF5/storm_d1_mass_<CASE>.h5', '/d1_mass' }
    { 'FILTERS/all_500_TSD_SAL_DUST.h5' '/filter' 'HDF5/<CASE>/HDF5/d2_mass-<CASE>-AS-2006-08-20-120000-g3.h5' '/dust2_mass' 'VAPOR/HDF5/storm_d2_mass_<CASE>.h5', '/d2_mass' }

    { 'FILTERS/all_500_TSD_SAL_DUST.h5' '/filter' 'HDF5/<CASE>/HDF5/d1_num-<CASE>-AS-2006-08-20-120000-g3.h5' '/dust1_concen' 'VAPOR/HDF5/storm_d1_num_<CASE>.h5', '/d1_num' }
    { 'FILTERS/all_500_TSD_SAL_DUST.h5' '/filter' 'HDF5/<CASE>/HDF5/d2_num-<CASE>-AS-2006-08-20-120000-g3.h5' '/dust2_concen' 'VAPOR/HDF5/storm_d2_num_<CASE>.h5', '/d2_num' }
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
      h5create(OutFile, Ovname, Vsize, 'ChunkSize', Csize, 'Deflate', 6, 'Shuffle', true);

      % prep for reading input files one timestep at a time
      FILTER_DS = ncgeodataset(FilterFile);
      DUST_DS   = ncgeodataset(DustFile);

      FILTER_VAR = FILTER_DS.geovariable(Fvname);
      DUST_VAR   = DUST_DS.geovariable(Dvname);

      % Do one time step at a time to help keep memory requirements lower.
      %
      % DUST_VAR has size of (Nx,Ny,Nz,Nt) and FILTER_VAR has size of (Nx,Ny,1,Nt) so
      % the single in FILTER_VAR needs to be replicated Nz times.
      for it = 1:Nt
        % nctoolbox, for a single timestep, returns dimensions as
        % (z,y,x) or (y,x). Run permute to put these into (x,y,z) and (x,y) order.
        FILTER = squeeze(FILTER_VAR.data(it,:,:,:));
        DUST   = squeeze(DUST_VAR.data(it,:,:,:));

        FILTER = permute(FILTER, [ 2 1 ]);
        DUST = permute(DUST, [ 3 2 1 ]);

        % Replicate levels
        FILTER = repmat(FILTER, [ 1 1 Nz ]);

        % Select dust only where filter has 1's (filter entries are either 1's or 0's)
        DUST = DUST .* FILTER;

        % Write DUST to output file. 
        Dstart = [ 1 1 1 it ];
        Dcount = [ Nx Ny Nz 1 ];
        h5write (OutFile, Ovname, DUST, Dstart, Dcount);

        if (mod(it, 10) == 0)
          fprintf('  Completed timestep: %d\n', it);
        end
      end
      
      % Attach dimensions
      AttachDimensionsXyzt(OutFile, Ovname, DimOrder, Xname, Yname, Zname, Tname);
    end
  end
end 
