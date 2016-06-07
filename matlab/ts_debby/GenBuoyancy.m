function [ ] = GenBuoyancy()
% GenBuoyancy function to calculate buoyancy from theta rho values

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
    { 'FILTERS/core_all_<CASE>.h5' '/filter' 'HDF5/<CASE>/HDF5/theta_rho-<CASE>-AS-2006-08-20-120000-g3.h5' '/theta_rho' 'HDF5/<CASE>/HDF5/core_buoy_acc-<CASE>-AS-2006-08-20-120000-g3.h5' '/buoy_acc' 'm/s2' 'buoyancy_acceleration' }
    { 'FILTERS/rb_all_<CASE>.h5' '/filter' 'HDF5/<CASE>/HDF5/theta_rho-<CASE>-AS-2006-08-20-120000-g3.h5' '/theta_rho' 'HDF5/<CASE>/HDF5/rband_buoy_acc-<CASE>-AS-2006-08-20-120000-g3.h5' '/buoy_acc' 'm/s2' 'buoyancy_acceleration' }
    };
  Nfiles = length(FileList);


  % For output notations
  DimOrder = { 'x' 'y' 'z' 't' };
  DimNames = 'x y z t';
  Units    = 'm/s2';
  LongName = 'buoyancy_accelleration';

  for icase = 1:Ncases
    Case = CaseList{icase};

    fprintf('*****************************************************************\n');
    fprintf('Generating  anomalies:\n');
    fprintf('  Case: %s\n', Case);
    fprintf('\n');

    for ifile = 1:Nfiles
      Ffile    = FileList{ifile}{1};
      Fvname   = FileList{ifile}{2};
      Tfile    = FileList{ifile}{3};
      Tvname   = FileList{ifile}{4};
      Ofile    = FileList{ifile}{5};
      Ovname   = FileList{ifile}{6};

      FilterFile   = regexprep(Ffile, '<CASE>', Case);
      ThetaRhoFile = regexprep(Tfile, '<CASE>', Case);
      OutFile      = regexprep(Ofile, '<CASE>', Case);

      fprintf('  Reading: %s (%s)\n', FilterFile,   Fvname);
      fprintf('  Reading: %s (%s)\n', ThetaRhoFile, Tvname);
      fprintf('\n');
      fprintf('  Writing: %s (%s)\n', OutFile, Ovname);
      fprintf('\n');

      % Read in coordinates, set size of input data
      X = squeeze(h5read(ThetaRhoFile, '/x_coords'));
      Y = squeeze(h5read(ThetaRhoFile, '/y_coords'));
      Z = squeeze(h5read(ThetaRhoFile, '/z_coords'));
      T = squeeze(h5read(ThetaRhoFile, '/t_coords'));

      Nx = length(X);
      Ny = length(Y);
      Nz = length(Z);
      Nt = length(T);

      % get old file out of the way
      if (exist(OutFile, 'file') == 2)
        delete(OutFile);
      end

      % Create the output dataset, give time dimension ability to continuously expand
      Vsize = [ Nx Ny Nz Inf ];
      h5create(OutFile, Ovname, Vsize, 'DataType', 'single', 'ChunkSize', [ Nx Ny 1 1 ], 'Deflate', 1, 'Shuffle', 1);

      % Do one time step at a time in order to not overwhelm the memory usage.
      for it = 1:Nt
        Start = [ 1 1 1 it ];
        Fcount = [ Nx Ny  1 1 ];
        Tcount = [ Nx Ny Nz 1 ];
        FDATA = squeeze(h5read(FilterFile, Fvname, Start, Fcount));
        TDATA = squeeze(h5read(ThetaRhoFile, Tvname, Start, Tcount));

        % Repeat FDATA along the z dimension so that it's dimensions match those of TDATA
        FDATA = repmat(FDATA, [ 1 1 Nz ]);

        % Apply selection using the filter data
        TH_RHO = FDATA .* TDATA;
        TH_RHO(TH_RHO == 0) = nan;

        % Find the mean theta rho, ignoring places where nans (originally zero values) exist
        % Horizontal mean --> apply nanmean over first two dimensions
        % Then repeat the matrix out to (x,y,z) for subsequent calculations
        TH_RHO_BAR = nanmean(TH_RHO, 1);
        TH_RHO_BAR = nanmean(TH_RHO_BAR, 2);
        TH_RHO_BAR = repmat(TH_RHO_BAR, [ Nx Ny 1 ]);

        % Acceleration due to buoyancy is calculated:
        %
        %   B  = g (TH_RHO - TH_RHO_BAR) / (TH_RHO_BAR);
        %
        %       g = 9.8 m/s^2
        BUOY = 9.8 .* (TH_RHO - TH_RHO_BAR) ./ TH_RHO_BAR;

        % Write the data for this time step
        h5write (OutFile, Ovname, BUOY, Start, Tcount);

        % report progress to user
        if (mod(it,10) == 0)
          fprintf('    Working, time step = %d\n', it);
        end
      end
      fprintf('\n');

      % Write coordinates into the output file
      Xname = '/x_coords';
      Yname = '/y_coords';
      Zname = '/z_coords';
      Tname = '/t_coords';

      CreateDimensionsXyzt(OutFile, X, Y, Z, T, Xname, Yname, Zname, Tname);

      % Add COARDS annotations on coordinates
      NotateDimensionsXyzt(OutFile, Xname, Yname, Zname, Tname);

      % Attach dimensions
      AttachDimensionsXyzt(OutFile, Ovname, DimOrder, Xname, Yname, Zname, Tname);

      % Put in COARDS annotations on variable
      NotateVariableXyzt(OutFile, Ovname, Units, LongName, DimNames);

    end
  end
end 
