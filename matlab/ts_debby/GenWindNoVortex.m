function [ ] = GenWindNoVortex()
% GenWindNoVortex - generate u, v winds that exclude the TS Debby vortex
%
%  This script will remove the vortex in the u,v winds by:
%    1) use the all_500 filter to identify vortex extent
%    2) remove data within vortex extent
%    3) use inpaintn to interpolate the removed data from
%       the surrounding environmental winds
%

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
    fprintf('Generating wind data, removing vortex: %s\n', Case);
    fprintf('\n');

    UinFile = sprintf('HDF5/%s/HDF5/u_lite-%s-AS-2006-08-20-120000-g3.h5', Case, Case);
    VinFile = sprintf('HDF5/%s/HDF5/v_lite-%s-AS-2006-08-20-120000-g3.h5', Case, Case);

    FinFile = sprintf('FILTERS/all_500_lite_%s.h5', Case);

    Uvname = '/u';
    Vvname = '/v';
    Fvname = '/filter';

    fprintf('  Reading: %s (%s)\n', UinFile, Uvname);
    fprintf('  Reading: %s (%s)\n', VinFile, Vvname);
    fprintf('  Reading: %s (%s)\n', FinFile, Fvname);
    fprintf('\n');

    % Process one time step at a time
    X = squeeze(h5read(UinFile, Xvname));
    Y = squeeze(h5read(UinFile, Yvname));
    Z = squeeze(h5read(UinFile, Zvname));
    T = squeeze(h5read(UinFile, Tvname));

    Nx = length(X);
    Ny = length(Y);
    Nz = length(Z);
    Nt = length(T);
 
    % Create the output files and dataset so that output can be written one time step at a time.
    UoutFile = sprintf('HDF5/%s/HDF5/no_vortex_u_lite-%s-AS-2006-08-20-120000-g3.h5', Case, Case);
    VoutFile = sprintf('HDF5/%s/HDF5/no_vortex_v_lite-%s-AS-2006-08-20-120000-g3.h5', Case, Case);

    if (exist(UoutFile, 'file') == 2)
      delete(UoutFile);
    end
    if (exist(VoutFile, 'file') == 2)
      delete(VoutFile);
    end

    fprintf('  Writing %s (%s)\n', UoutFile, Uvname);
    fprintf('  Writing %s (%s)\n', VoutFile, Vvname);
    fprintf('\n');

    % copy the dimension vars into the output file
    CreateDimensionsXyzt(UoutFile, X, Y, Z, T, Xvname, Yvname, Zvname, Tvname);
    CreateDimensionsXyzt(VoutFile, X, Y, Z, T, Xvname, Yvname, Zvname, Tvname);

    NotateDimensionsXyzt(UoutFile, Xvname, Yvname, Zvname, Tvname);
    NotateDimensionsXyzt(VoutFile, Xvname, Yvname, Zvname, Tvname);

    % create output datasets that can be written one time step at a time
    Vsize = [ Nx Ny Nz Inf ];
    Csize = [ Nx Ny Nz   1 ];
    h5create(UoutFile, Uvname, Vsize, 'DataType', 'single', 'ChunkSize', Csize, 'Deflate', 6, 'Shuffle', true);
    h5create(VoutFile, Vvname, Vsize, 'DataType', 'single', 'ChunkSize', Csize, 'Deflate', 6, 'Shuffle', true);

    for it = 1:Nt
      Start     = [ 1 1 1 it ];
      IoCount   = [ Nx Ny Nz 1 ];
      FiltCount = [ Nx Ny 1 1 ];

      U = squeeze(h5read(UinFile, Uvname, Start, IoCount));
      V = squeeze(h5read(VinFile, Vvname, Start, IoCount));

      FILTER = squeeze(h5read(FinFile, Fvname, Start, FiltCount));
      FILTER_3D = repmat(FILTER, [ 1 1 Nz ]);

      % Replace entries in U and V with nans wherever filter is a '1'
      U(FILTER_3D == 1) = nan;
      V(FILTER_3D == 1) = nan;

      % Interpolate level by level
      U_NV = zeros([ Nx Ny Nz ]);
      V_NV = zeros([ Nx Ny Nz ]);
      for k = 1:Nz
        U_NV(:,:,k) = single(inpaintn(squeeze(U(:,:,k))));
        V_NV(:,:,k) = single(inpaintn(squeeze(V(:,:,k))));
      end

      h5write(UoutFile, Uvname, U_NV, Start, IoCount);
      h5write(VoutFile, Vvname, V_NV, Start, IoCount);

      if (mod(it, 10) == 0)
        fprintf('  Completed timestep: %d\n', it);
      end
    end
    fprintf('  Processed %d timesteps\n', it);
    fprintf('\n');

    % attach dimensions to output u,v for GRADS
    DimOrder = { 'x' 'y' 'z' 't' };
    AttachDimensionsXyzt(UoutFile, Uvname, DimOrder, Xvname, Yvname, Zvname, Tvname);
    AttachDimensionsXyzt(VoutFile, Vvname, DimOrder, Xvname, Yvname, Zvname, Tvname);

  end % cases
end % function