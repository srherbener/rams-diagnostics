function [ ] = GenWindShearMeas()
% GenWindShearMeas - create time series of average shear magnitude from wind shear field

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

  for icase = 1:Ncases
    Case = CaseList{icase};

    fprintf('*****************************************************************\n');
    fprintf('Generating average wind shear: %s\n', Case);
    fprintf('  No vortex: %d\n', NoVortex);
    fprintf('\n');

    if (NoVortex == 1)
      InFile = sprintf('HDF5/%s/HDF5/shear_nv_lite-%s-AS-2006-08-20-120000-g3.h5', Case, Case);
    else
      InFile = sprintf('HDF5/%s/HDF5/shear_lite-%s-AS-2006-08-20-120000-g3.h5', Case, Case);
    end
    FiltFile = sprintf('FILTERS/all_500_lite_%s.h5', Case);

    InVname = '/mag_shear';
    FiltVname = '/filter';

    fprintf('  Reading: %s (%s)\n', InFile, InVname);
    fprintf('  Reading: %s (%s)\n', FiltFile, FiltVname);
    fprintf('\n');

    % Process one time step at a time
    % Find the average u and average v at 850mb and 200mb levels
    % Then take vector difference to form the shear vector
    
    X = squeeze(h5read(InFile, Xvname));
    Y = squeeze(h5read(InFile, Yvname));
    Z = squeeze(h5read(InFile, Zvname));
    T = squeeze(h5read(InFile, Tvname));

    Nx = length(X);
    Ny = length(Y);
    Nz = length(Z);
    Nt = length(T);

    MAG_SHEAR_AVG = zeros([ Nt 1 ]);

    for it = 1:Nt
      Start = [  1  1 it ];
      Count = [ Nx Ny  1 ];
      MAG_SHEAR = squeeze(h5read(InFile, InVname, Start, Count));

      Start = [  1  1  1 it ];
      Count = [ Nx Ny  1  1 ];
      FILTER    = squeeze(h5read(FiltFile, FiltVname, Start, Count));

      % Only want the average performed over the 500 km radius defined
      % in FILTER. FILTER has 1's within the 500 radius and 0's elsewhere,
      % therefore can use FILTER to select out of MAG_SHEAR.
      MAG_SHEAR_AVG(it) = mean(MAG_SHEAR(FILTER == 1));

      if (mod(it, 10) == 0)
        fprintf('  Completed timestep: %d\n', it);
      end
    end
    fprintf('  Processed %d timesteps\n', it);
    fprintf('\n');

    % Output
    if (NoVortex == 1)
      OutFile = sprintf('DIAGS/wind_shear_avg_no_vortex_%s.h5', Case);
    else
      OutFile = sprintf('DIAGS/wind_shear_avg_%s.h5', Case);
    end
    if (exist(OutFile, 'file') == 2)
      delete(OutFile);
    end

    CreateDimensionsXyzt(OutFile, X, Y, Z, T, Xvname, Yvname, Zvname, Tvname);
    NotateDimensionsXyzt(OutFile, Xvname, Yvname, Zvname, Tvname);

    Vsize = Nt;
    DimOrder = { 't' };

    Ovname = '/mag_shear';
    fprintf('  Writing %s (%s)\n', OutFile, Ovname);
    h5create(OutFile, Ovname, Vsize);
    h5write(OutFile, Ovname, MAG_SHEAR_AVG);
    AttachDimensionsXyzt(OutFile, Ovname, DimOrder, Xvname, Yvname, Zvname, Tvname);

  end % cases
end % function
