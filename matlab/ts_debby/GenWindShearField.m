function [ ] = GenWindShearField()
% GenWindShearField - generate wind shear vector field
%
% Use method as described in Dunion and Velden, 2004
%   average wind in 150-350mb layer minus the average wind in 700-925mb layer
%
% Using dry calculations (adiabatic lapse rate), conversion of pressure to height yields:
%    150 mb   13,600 m
%    350 mb    8,120 m
%    700 mb    3,010 m
%    925 mb      762 m
%
% Create a field by taking average u,v along columns in the upper and lower layers, and
% then calculating shear at each grid point.
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

  % NoVortex
  %   1 -> use the u,v data with the vortex subtracted out
  %   0 -> use the original u,v data
  NoVortex = 1;

  % approximate heights for the measurement layers
  Z150 = 13600; % m
  Z350 =  8120; % m
  Z700 =  3010; % m
  Z925 =   762; % m

  for icase = 1:Ncases
    Case = CaseList{icase};

    fprintf('*****************************************************************\n');
    fprintf('Generating wind shear field: %s\n', Case);
    fprintf('  No vortex: %d\n', NoVortex);
    fprintf('\n');

    if (NoVortex == 1)
      Uvname = '/u_no_vortex';
      Vvname = '/v_no_vortex';
      Ufile = sprintf('HDF5/%s/HDF5/no_vortex_lite-%s-AS-2006-08-20-120000-g3.h5', Case, Case);
      Vfile = sprintf('HDF5/%s/HDF5/no_vortex_lite-%s-AS-2006-08-20-120000-g3.h5', Case, Case);
    else
      Uvname = '/u';
      Vvname = '/v';
      Ufile = sprintf('HDF5/%s/HDF5/u_lite-%s-AS-2006-08-20-120000-g3.h5', Case, Case);
      Vfile = sprintf('HDF5/%s/HDF5/v_lite-%s-AS-2006-08-20-120000-g3.h5', Case, Case);
    end

    fprintf('  Reading: %s (%s)\n', Ufile, Uvname);
    fprintf('  Reading: %s (%s)\n', Vfile, Vvname);
    fprintf('\n');

    % Process one time step at a time
    % Create mean column u and v for the two layers
    % Then take vector difference to form the shear vector
    
    X = squeeze(h5read(Ufile, Xvname));
    Y = squeeze(h5read(Ufile, Yvname));
    Z = squeeze(h5read(Ufile, Zvname));
    T = squeeze(h5read(Ufile, Tvname));

    Nx = length(X);
    Ny = length(Y);
    Nz = length(Z);
    Nt = length(T);

    % Find the indices for the measurement layers
    ZL1 = find(Z >= Z925, 1, 'first');
    ZL2 = find(Z <= Z700, 1, 'last');
    ZU1 = find(Z >= Z350, 1, 'first');
    ZU2 = find(Z <= Z150, 1, 'last');

    Nzlow = (ZL2 - ZL1) + 1;
    Nzhigh = (ZU2 - ZU1) + 1;

    % Output file, set up to be able to write output one time step at a time
    % Output will not have z-dimension, ie (x,y,t)
    if (NoVortex == 1)
      OutFile = sprintf('HDF5/%s/HDF5/shear_nv_lite-%s-AS-2006-08-20-120000-g3.h5', Case, Case);
    else
      OutFile = sprintf('HDF5/%s/HDF5/shear_lite-%s-AS-2006-08-20-120000-g3.h5', Case, Case);
    end
    if (exist(OutFile, 'file') == 2)
      delete(OutFile);
    end

    CreateDimensionsXyzt(OutFile, X, Y, Z, T, Xvname, Yvname, Zvname, Tvname);
    NotateDimensionsXyzt(OutFile, Xvname, Yvname, Zvname, Tvname);

    Vsize = [ Nx Ny Inf ];
    Csize = [ Nx Ny   1 ];

    USvname = '/u_shear';
    fprintf('  Writing %s (%s)\n', OutFile, USvname);
    h5create(OutFile, USvname, Vsize, 'DataType', 'single', 'ChunkSize', Csize, 'Deflate', 6, 'Shuffle', true);

    VSvname = '/v_shear';
    fprintf('  Writing %s (%s)\n', OutFile, VSvname);
    h5create(OutFile, VSvname, Vsize, 'DataType', 'single', 'ChunkSize', Csize, 'Deflate', 6, 'Shuffle', true);

    MSvname = '/mag_shear';
    fprintf('  Writing %s (%s)\n', OutFile, MSvname);
    h5create(OutFile, MSvname, Vsize, 'DataType', 'single', 'ChunkSize', Csize, 'Deflate', 6, 'Shuffle', true);

    ASvname = '/ang_shear';
    fprintf('  Writing %s (%s)\n', OutFile, ASvname);
    h5create(OutFile, ASvname, Vsize, 'DataType', 'single', 'ChunkSize', Csize, 'Deflate', 6, 'Shuffle', true);

    fprintf('\n');

    % work on time step at a time in order to reduce memory requirements
    for it = 1:Nt
      Start = [ 1 1 ZL1 it ];
      Count = [ Nx Ny Nzlow 1 ];
      ULOW = squeeze(h5read(Ufile, Uvname, Start, Count));
      VLOW = squeeze(h5read(Vfile, Vvname, Start, Count));

      Start = [ 1 1 ZU1 it ];
      Count = [ Nx Ny Nzhigh 1 ];
      UHIGH = squeeze(h5read(Ufile, Uvname, Start, Count));
      VHIGH = squeeze(h5read(Vfile, Vvname, Start, Count));

      % Form the column average (z dimension is number 3)
      ULOW = squeeze(mean(ULOW, 3));
      VLOW = squeeze(mean(VLOW, 3));
      UHIGH = squeeze(mean(UHIGH, 3));
      VHIGH = squeeze(mean(VHIGH, 3));

      % shear - cartesian
      U_SHEAR = UHIGH - ULOW;
      V_SHEAR = VHIGH - VLOW;

      % shear - mag/angle
      MAG_SHEAR = sqrt(U_SHEAR.^2 + V_SHEAR.^2);
      ANG_SHEAR = atan2(V_SHEAR, U_SHEAR);

      % write out this timestep
      Start = [  1  1 it ];
      Count = [ Nx Ny  1 ];
      h5write(OutFile, USvname, U_SHEAR, Start, Count);
      h5write(OutFile, VSvname, V_SHEAR, Start, Count);
      h5write(OutFile, MSvname, MAG_SHEAR, Start, Count);
      h5write(OutFile, ASvname, ANG_SHEAR, Start, Count);

      if (mod(it, 10) == 0)
        fprintf('    Completed timestep: %d\n', it);
      end
    end
    fprintf('    Processed %d timesteps\n', it);
    fprintf('\n');

    % attach dimensions (for GRADS)
    DimOrder = { 'x' 'y' 't' };
    AttachDimensionsXyzt(OutFile, USvname, DimOrder, Xvname, Yvname, Zvname, Tvname);
    AttachDimensionsXyzt(OutFile, VSvname, DimOrder, Xvname, Yvname, Zvname, Tvname);
    AttachDimensionsXyzt(OutFile, MSvname, DimOrder, Xvname, Yvname, Zvname, Tvname);
    AttachDimensionsXyzt(OutFile, ASvname, DimOrder, Xvname, Yvname, Zvname, Tvname);

  end % cases
end % function
