function [ ] = GenMoistRamsInit()
% GenMoistRamsInit function to generate moist environment data for TS Debby RAMS simualations

  % Read in RAMS moist sim info (created with GenMoistRamsInfo) and select out the
  % regions of RH that the masks indicate.

  Ngrids = 3;
  Hdir = 'HDF5';

  InFile  = sprintf('%s/MoistRamsInfo.h5', Hdir);
  OutFile = sprintf('%s/TS_Debby_Init.h5', Hdir);

  fprintf('*********************************************************************\n');
  fprintf('Selecting new RH values for RAMS TS Debby moist simulation initialization:\n');
  fprintf('\n');

  fprintf('  Writing init data into: %s\n', OutFile);
  fprintf('\n');

  hdf5write(OutFile, 'Header', 'TS Debby Moist Initialization');

  % Read in each grid and interpolate the mask data to that grid. Also, calculate the new RH
  % values from the data in the var files. Write out the interpolated masks and RH values for
  % processing into initialization data for RAMS.

  for igrid = 1:Ngrids
    MaskVar = sprintf('RAMS_MASK_G%d', igrid);
    RhVar   = sprintf('RAMS_RH_G%d',   igrid);

    fprintf('  Grid: %d\n', igrid);
    fprintf('    Reading RAMS info file: %s -> %s, %s\n', InFile, MaskVar, RhVar);
    fprintf('\n');

    MASK = squeeze(hdf5read(InFile, MaskVar));
    RH   = squeeze(hdf5read(InFile, RhVar));

    [ Nx Ny Nz ] = size(RH);

    % Create an (i,j,k) grid mesh. Then select, using MASK, out of RH and
    % the I, J, K vectors to create lists of corresponding i,j,k,rh values.
    [ I J K ] = ndgrid(1:Nx, 1:Ny, 1:Nz);

    % This will create linear arrays (vectors) which is ideal since the RAMS hdf5
    % input routine (shdf5_irec) can be used. This also minimizes the amount of data
    % since only the selected points are saved.
    I_VALS  = I(MASK == 1);
    J_VALS  = J(MASK == 1);
    K_VALS  = K(MASK == 1);
    RH_VALS = RH(MASK == 1);

    % Don't allow RH to go above 1. Ie, don't allow supersaturation since this could artificially
    % create clouds.
    RH_VALS(RH_VALS > 0.999) = 0.999;

    % Output separate vectors for each grid
    OutVar = sprintf('G%d_I', igrid);
    hdf5write(OutFile, OutVar, I_VALS, 'WriteMode', 'append');
    OutVar = sprintf('G%d_J', igrid);
    hdf5write(OutFile, OutVar, J_VALS, 'WriteMode', 'append');
    OutVar = sprintf('G%d_K', igrid);
    hdf5write(OutFile, OutVar, K_VALS, 'WriteMode', 'append');
    OutVar = sprintf('G%d_RH', igrid);
    hdf5write(OutFile, OutVar, RH_VALS, 'WriteMode', 'append');
  end
end
