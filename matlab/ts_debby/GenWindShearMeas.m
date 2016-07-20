function [ ] = GenWindShearMeas()
% GenWindShearMeas - generate wind shear vector, magnitude, angle in context of SAL
%
% Use method as described in Dunion and Velden, 2004
%   average wind in 150-350mb layer minus the average wind in 700-925mb layer
%
% Using dry calculations (adiabatic lapse rate), conversion of pressure to height yields:
%    150 mb   13,600 m
%    350 mb    8,120 m
%    700 mb    3,010 m
%    925 mb      76s m
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
  NoVortex = 0;

  % approximate heights for the measurement layers
  Z150 = 13600; % m
  Z350 =  8120; % m
  Z700 =  3010; % m
  Z925 =   762; % m

  for icase = 1:Ncases
    Case = CaseList{icase};

    fprintf('*****************************************************************\n');
    fprintf('Generating wind shear vector: %s\n', Case);
    fprintf('  No vortex: %d\n', NoVortex);
    fprintf('\n');

    if (NoVortex == 1)
      Ufile = sprintf('HDF5/%s/HDF5/no_vortex_u_lite-%s-AS-2006-08-20-120000-g3.h5', Case, Case);
      Vfile = sprintf('HDF5/%s/HDF5/no_vortex_v_lite-%s-AS-2006-08-20-120000-g3.h5', Case, Case);
    else
      Ufile = sprintf('HDF5/%s/HDF5/u_lite-%s-AS-2006-08-20-120000-g3.h5', Case, Case);
      Vfile = sprintf('HDF5/%s/HDF5/v_lite-%s-AS-2006-08-20-120000-g3.h5', Case, Case);
    end
    Ffile = sprintf('FILTERS/all_500_lite_%s.h5', Case);

    Uvname = '/u';
    Vvname = '/v';
    Fvname = '/filter';

    fprintf('  Reading: %s (%s)\n', Ufile, Uvname);
    fprintf('  Reading: %s (%s)\n', Vfile, Vvname);
    fprintf('  Reading: %s (%s)\n', Ffile, Fvname);
    fprintf('\n');

    % Process one time step at a time
    % Find the average u and average v at 850mb and 200mb levels
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

    ULOW_AVG = zeros([ Nt 1 ]);
    VLOW_AVG = zeros([ Nt 1 ]);
    UHIGH_AVG = zeros([ Nt 1 ]);
    VHIGH_AVG = zeros([ Nt 1 ]);

    U_SHEAR = zeros([ Nt 1 ]);
    V_SHEAR = zeros([ Nt 1 ]);

    MAG_SHEAR = zeros([ Nt 1 ]);
    ANGLE_SHEAR = zeros([ Nt 1 ]);

    for it = 1:Nt
      Start = [ 1 1 ZL1 it ];
      Count = [ Nx Ny Nzlow 1 ];
      ULOW = squeeze(h5read(Ufile, Uvname, Start, Count));
      VLOW = squeeze(h5read(Vfile, Vvname, Start, Count));

      Start = [ 1 1 ZU1 it ];
      Count = [ Nx Ny Nzhigh 1 ];
      UHIGH = squeeze(h5read(Ufile, Uvname, Start, Count));
      VHIGH = squeeze(h5read(Vfile, Vvname, Start, Count));

      Start = [ 1 1 1 it ];
      Count = [ Nx Ny 1 1 ];
      FILTER = squeeze(h5read(Ffile, Fvname, Start, Count));

      FILTER_LOW = repmat(FILTER, [ 1 1 Nzlow ]);
      FILTER_HIGH = repmat(FILTER, [ 1 1 Nzhigh ]);

      % Only want the average performed over the 500 km radius defined
      % in FILTER. FILTER has 1's within the 500 radius and 0's elsewhere,
      % therefore:
      %    multiply u, v by filter
      %    change zeros to nans
      %    use nanmean to do the average
      ULOW = ULOW .* FILTER_LOW;
      ULOW(ULOW == 0) = nan;

      VLOW = VLOW .* FILTER_LOW;
      VLOW(VLOW == 0) = nan;

      UHIGH = UHIGH .* FILTER_HIGH;
      UHIGH(UHIGH == 0) = nan;

      VHIGH = VHIGH .* FILTER_HIGH;
      VHIGH(VHIGH == 0) = nan;

      % write out average u, v values into the output file to
      % help inspect results 
      ULOW_AVG(it) = nanmean(ULOW(:));
      VLOW_AVG(it) = nanmean(VLOW(:));
      UHIGH_AVG(it) = nanmean(UHIGH(:));
      VHIGH_AVG(it) = nanmean(VHIGH(:));

      U_SHEAR(it) = UHIGH_AVG(it) - ULOW_AVG(it);
      V_SHEAR(it) = VHIGH_AVG(it) - VLOW_AVG(it);

      if (mod(it, 10) == 0)
        fprintf('  Completed timestep: %d\n', it);
      end
    end
    fprintf('  Processed %d timesteps\n', it);
    fprintf('\n');

    MAG_SHEAR = sqrt(U_SHEAR.^2 + V_SHEAR.^2);
    ANGLE_SHEAR = atan2(V_SHEAR, U_SHEAR);

    % Output
    if (NoVortex == 1)
      OutFile = sprintf('DIAGS/wind_shear_no_vortex_%s.h5', Case);
    else
      OutFile = sprintf('DIAGS/wind_shear_%s.h5', Case);
    end
    if (exist(OutFile, 'file') == 2)
      delete(OutFile);
    end

    CreateDimensionsXyzt(OutFile, X, Y, Z, T, Xvname, Yvname, Zvname, Tvname);
    NotateDimensionsXyzt(OutFile, Xvname, Yvname, Zvname, Tvname);

    Vsize = Nt;
    DimOrder = { 't' };

    Ovname = '/u_low_avg';
    fprintf('  Writing %s (%s)\n', OutFile, Ovname);
    h5create(OutFile, Ovname, Vsize);
    h5write(OutFile, Ovname, ULOW_AVG);
    AttachDimensionsXyzt(OutFile, Ovname, DimOrder, Xvname, Yvname, Zvname, Tvname);

    Ovname = '/v_low_avg';
    fprintf('  Writing %s (%s)\n', OutFile, Ovname);
    h5create(OutFile, Ovname, Vsize);
    h5write(OutFile, Ovname, VLOW_AVG);
    AttachDimensionsXyzt(OutFile, Ovname, DimOrder, Xvname, Yvname, Zvname, Tvname);

    Ovname = '/u_high_avg';
    fprintf('  Writing %s (%s)\n', OutFile, Ovname);
    h5create(OutFile, Ovname, Vsize);
    h5write(OutFile, Ovname, UHIGH_AVG);
    AttachDimensionsXyzt(OutFile, Ovname, DimOrder, Xvname, Yvname, Zvname, Tvname);

    Ovname = '/v_high_avg';
    fprintf('  Writing %s (%s)\n', OutFile, Ovname);
    h5create(OutFile, Ovname, Vsize);
    h5write(OutFile, Ovname, VHIGH_AVG);
    AttachDimensionsXyzt(OutFile, Ovname, DimOrder, Xvname, Yvname, Zvname, Tvname);

    Ovname = '/u_shear';
    fprintf('  Writing %s (%s)\n', OutFile, Ovname);
    h5create(OutFile, Ovname, Vsize);
    h5write(OutFile, Ovname, U_SHEAR);
    AttachDimensionsXyzt(OutFile, Ovname, DimOrder, Xvname, Yvname, Zvname, Tvname);

    Ovname = '/v_shear';
    fprintf('  Writing %s (%s)\n', OutFile, Ovname);
    h5create(OutFile, Ovname, Vsize);
    h5write(OutFile, Ovname, V_SHEAR);
    AttachDimensionsXyzt(OutFile, Ovname, DimOrder, Xvname, Yvname, Zvname, Tvname);

    Ovname = '/mag_shear';
    fprintf('  Writing %s (%s)\n', OutFile, Ovname);
    h5create(OutFile, Ovname, Vsize);
    h5write(OutFile, Ovname, MAG_SHEAR);
    AttachDimensionsXyzt(OutFile, Ovname, DimOrder, Xvname, Yvname, Zvname, Tvname);

    Ovname = '/angle_shear';
    fprintf('  Writing %s (%s)\n', OutFile, Ovname);
    h5create(OutFile, Ovname, Vsize);
    h5write(OutFile, Ovname, ANGLE_SHEAR);
    AttachDimensionsXyzt(OutFile, Ovname, DimOrder, Xvname, Yvname, Zvname, Tvname);

  end % cases
end % function
