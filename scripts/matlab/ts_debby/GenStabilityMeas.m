function [ ] = GenStabilityMeas()
% GenStabilityMeas - generate atmospheric stability measurements
%
% This script will generate Brunt-Vaisala Frequency (N^2) values
% from input virtual potential temperature:
%
%   N^2 = (g/theta_v) * (dtheta_v/dz)
%

  g = 9.81; % gravity acceleration, m/s

  Cp = 1005; % heat capacity of air at constant pressure, J/(kg K)
  Lv = 2.5e6; % latent heat of vaporization for water, J/kg

  % list of simulation cases
  CaseList = {
    'TSD_SAL_DUST'
    'TSD_SAL_NODUST'
    'TSD_NONSAL_DUST'
    'TSD_NONSAL_NODUST'
    };
  Ncases = length(CaseList);

  Xvname = '/x_coords';
  Yvname = '/y_coords';
  Zvname = '/z_coords';
  Tvname = '/t_coords';

  for icase = 1:Ncases
    Case = CaseList{icase};

    fprintf('*****************************************************************\n');
    fprintf('Generating stability measurements: %s\n', Case);
    fprintf('\n');

    ThvFile = sprintf('HDF5/%s/HDF5/theta_v_lite-%s-AS-2006-08-20-120000-g3.h5', Case, Case);
    ThvVname = '/theta_v';

    VapFile = sprintf('HDF5/%s/HDF5/vapor_lite-%s-AS-2006-08-20-120000-g3.h5', Case, Case);
    VapVname = '/vapor';

    TempFile = sprintf('HDF5/%s/HDF5/tempc_lite-%s-AS-2006-08-20-120000-g3.h5', Case, Case);
    TempVname = '/tempc';

    fprintf('  Reading: %s (%s)\n', ThvFile, ThvVname);
    fprintf('  Reading: %s (%s)\n', VapFile, VapVname);
    fprintf('  Reading: %s (%s)\n', TempFile, TempVname);
    fprintf('\n');

    % Process one time step at a time
    X = squeeze(h5read(ThvFile, Xvname));
    Y = squeeze(h5read(ThvFile, Yvname));
    Z = squeeze(h5read(ThvFile, Zvname));
    T = squeeze(h5read(ThvFile, Tvname));

    Nx = length(X);
    Ny = length(Y);
    Nz = length(Z);
    Nt = length(T);

    % Want the resulting measurement values to be located at the input z levels.
    % To accomplish this, form z values at the k+1/2 levels which are then
    % used to calculated delta z at the k levels.
    Z_KHALF = zeros([ 1 Nz+1 ]);
    Z_KHALF(2:end-1) = (Z(2:end) + Z(1:end-1)) .* 0.5;
    Z_KHALF(1)   = Z(1) - (Z_KHALF(2) - Z(1));
    Z_KHALF(end) = Z(end) + (Z(end) - Z_KHALF(end-1));

    % Z_KHALF now contains the height values at the k+1/2 levels
    % that surround all of the original k levels. Note Z_KHALF has Nz+1
    % elements.
    DZ = Z_KHALF(2:end) - Z_KHALF(1:end-1);
    DZ = repmat(reshape(DZ, [ 1 1 Nz ]), [ Nx Ny 1 ]);

    % Need heights in 3D grid
    Z3D = repmat(reshape(Z, [ 1 1 Nz ]), [ Nx Ny 1 ]);

    % Create the output file and datasets so that output can be written one time step at a time.
    BvfFile = sprintf('HDF5/%s/HDF5/bvf_lite-%s-AS-2006-08-20-120000-g3.h5', Case, Case);
    BvfVname = '/brunt_vaisala_freq';
    if (exist(BvfFile, 'file') == 2)
      delete(BvfFile);
    end
    MseFile = sprintf('HDF5/%s/HDF5/mse_lite-%s-AS-2006-08-20-120000-g3.h5', Case, Case);
    MseVname = '/moist_static_energy';
    if (exist(MseFile, 'file') == 2)
      delete(MseFile);
    end

    % copy the dimension vars into the output file
    CreateDimensionsXyzt(BvfFile, X, Y, Z, T, Xvname, Yvname, Zvname, Tvname);
    NotateDimensionsXyzt(BvfFile, Xvname, Yvname, Zvname, Tvname);
    CreateDimensionsXyzt(MseFile, X, Y, Z, T, Xvname, Yvname, Zvname, Tvname);
    NotateDimensionsXyzt(MseFile, Xvname, Yvname, Zvname, Tvname);

    % create output datasets that can be written one time step at a time
    Vsize = [ Nx Ny Nz Inf ];
    Csize = [ Nx Ny Nz   1 ];

    h5create(BvfFile, BvfVname, Vsize, 'DataType', 'single', 'ChunkSize', Csize, 'Deflate', 6, 'Shuffle', true);
    h5create(MseFile, MseVname, Vsize, 'DataType', 'single', 'ChunkSize', Csize, 'Deflate', 6, 'Shuffle', true);

    fprintf('  Calculating:\n');
    fprintf('    Brunt-Vaisala Frequency\n');
    fprintf('    Moist Static Energy\n');
    for it = 1:Nt
      Start     = [ 1 1 1 it ];
      Count   = [ Nx Ny Nz 1 ];

      %########################################################
      % Brunt-Vaisala Frequency
      %########################################################
      THV = squeeze(h5read(ThvFile, ThvVname, Start, Count));

      % Calculate the output at the k vertical locations. Create a delta THV
      % array that contains values at the k+1/2 vertical locations. Similar
      % to DZ calculation above.
      THV_KHALF = zeros([ Nx Ny Nz+1 ]);
      THV_KHALF(:,:,2:end-1) = (THV(:,:,2:end) + THV(:,:,1:end-1)) .* 0.5;
      THV_KHALF(:,:,1)   = THV(:,:,1) - (THV_KHALF(:,:,2) - THV(:,:,1));
      THV_KHALF(:,:,end) = THV(:,:,end) + (THV(:,:,end) - THV_KHALF(:,:,end-1));

      DTHV = THV_KHALF(:,:,2:end) - THV_KHALF(:,:,1:end-1);

      NSQ = (g ./ THV) .* (DTHV ./ DZ);

      h5write(BvfFile, BvfVname, single(NSQ), Start, Count);

      %########################################################
      % Moist Static Energy
      %########################################################

      % Input vapor is in g/kg
      %   std units are kg/kg so multiply by 1e-3
      VAP = squeeze(h5read(VapFile, VapVname, Start, Count)) .* 1e-3;

      % Input temp is in deg C
      %   std unis are K so add 273.15
      TEMP = squeeze(h5read(TempFile, TempVname, Start, Count)) + 273.15;

      MSE = (Cp .* TEMP) + (g .* Z3D) + (Lv .* VAP); % Result is J/kg
      MSE = MSE .* 1e-3; % convert to kJ/kg

      h5write(MseFile, MseVname, single(MSE), Start, Count);

      if (mod(it, 10) == 0)
        fprintf('    Completed timestep: %d\n', it);
      end
    end
    fprintf('    Processed %d timesteps\n', it);
    fprintf('\n');

    % attach dimensions to output for GRADS
    fprintf('  Writing %s\n', BvfFile);
    fprintf('    %s\n', BvfVname);
    fprintf('  Writing %s\n', MseFile);
    fprintf('    %s\n', MseVname);

    DimOrder = { 'x' 'y' 'z' 't' };
    DimNames = 't z y x';

    Units = '1/s2';
    LongName = 'BruntVaisalFreq';
    AttachDimensionsXyzt(BvfFile, BvfVname, DimOrder, Xvname, Yvname, Zvname, Tvname);
    NotateVariableXyzt(BvfFile, BvfVname, Units, LongName, DimNames);

    Units = 'J/kg';
    LongName = 'MoistStaticEnergy';
    AttachDimensionsXyzt(MseFile, MseVname, DimOrder, Xvname, Yvname, Zvname, Tvname);
    NotateVariableXyzt(MseFile, MseVname, Units, LongName, DimNames);

  end % cases
end % function
