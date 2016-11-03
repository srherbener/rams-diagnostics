function [ ] = GenWindNoVortex()
% GenWindNoVortex - generate u, v winds that exclude the TS Debby vortex
%
%  This script will remove the vortex in the u,v winds by:
%
%  Method == 1
%    1) use the all_500 filter to identify vortex extent
%    2) remove data within vortex extent
%    3) use inpaintn to interpolate the removed data from
%       the surrounding environmental winds
%
%  Method == 2
%    1) use the all_500 filter to identify vortex extent
%    2) subtract out storm motion from vortex region
%    3) form circulation winds
%        a) azimuthal average
%    4) subtract out circulation from vortex region
%
%  Method == 3
%    Vortex separation as described in Kurihara et al. (1995, 1993)
%    1) Separate the total wind field into a basic field and
%       disturbance field. Use the smoothing method in Kurihara
%       et al. (1993) to do this.
%    2) Remove the vortex using a circular filter that tapers to zero at
%       at the circle edge in order to prevent a discontinuity in the
%       non-hurricane disturbance field when the vortex piece is removed.
%    3) The separated field then becomes the basic (smoothed) field plus
%       the non-hurricane piece of the disturbance field.

  Method = 3;

  R0 = 5; % degrees lat/lon (~500 km), for method 3
  L = R0 / 5; % for method 3

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
    fprintf('  Method: %d\n', Method);
    fprintf('\n');

    UinFile = sprintf('HDF5/%s/HDF5/u_lite-%s-AS-2006-08-20-120000-g3.h5', Case, Case);
    VinFile = sprintf('HDF5/%s/HDF5/v_lite-%s-AS-2006-08-20-120000-g3.h5', Case, Case);
    ScFile  = sprintf('HDF5/%s/HDF5/storm_center_lite-%s-AS-2006-08-20-120000-g3.h5',Case,Case);

    FinFile = sprintf('FILTERS/all_500_lite_%s.h5', Case);

    Uvname = '/u';
    Vvname = '/v';
    Fvname = '/filter';
    SsxVname = '/storm_speed_x';
    SsyVname = '/storm_speed_y';
    SxlocVname = '/press_cent_xloc';
    SylocVname = '/press_cent_yloc';

    UoutVname = '/u_orig';
    VoutVname = '/v_orig';
    UnvOutVname = '/u_no_vortex';
    VnvOutVname = '/v_no_vortex';

    fprintf('  Reading: %s (%s)\n', UinFile, Uvname);
    fprintf('  Reading: %s (%s)\n', VinFile, Vvname);
    fprintf('  Reading: %s (%s)\n', FinFile, Fvname);
    if (Method == 2)
      fprintf('  Reading: %s (%s, %s)\n', ScFile, SsxVname, SsyVname);
    end
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

    % Read in the storm speed and location components
    switch (Method)
      case 2
        StormSpeedX = squeeze(h5read(ScFile, SsxVname));
        StormSpeedY = squeeze(h5read(ScFile, SsyVname));

        StormLocX = squeeze(h5read(ScFile, SxlocVname));
        StormLocY = squeeze(h5read(ScFile, SylocVname));
      case 3
        StormLocX = squeeze(h5read(ScFile, SxlocVname));
        StormLocY = squeeze(h5read(ScFile, SylocVname));
    end
 
    % Create the output file and datasets so that output can be written one time step at a time.
    OutFile = sprintf('HDF5/%s/HDF5/no_vortex_lite-%s-AS-2006-08-20-120000-g3.h5', Case, Case);
    if (exist(OutFile, 'file') == 2)
      delete(OutFile);
    end

    % copy the dimension vars into the output file
    CreateDimensionsXyzt(OutFile, X, Y, Z, T, Xvname, Yvname, Zvname, Tvname);
    NotateDimensionsXyzt(OutFile, Xvname, Yvname, Zvname, Tvname);

    % create output datasets that can be written one time step at a time
    Vsize = [ Nx Ny Nz Inf ];
    Csize = [ Nx Ny Nz   1 ];

    h5create(OutFile, UoutVname, Vsize, 'DataType', 'single', 'ChunkSize', Csize, 'Deflate', 6, 'Shuffle', true);
    h5create(OutFile, VoutVname, Vsize, 'DataType', 'single', 'ChunkSize', Csize, 'Deflate', 6, 'Shuffle', true);
    h5create(OutFile, UnvOutVname, Vsize, 'DataType', 'single', 'ChunkSize', Csize, 'Deflate', 6, 'Shuffle', true);
    h5create(OutFile, VnvOutVname, Vsize, 'DataType', 'single', 'ChunkSize', Csize, 'Deflate', 6, 'Shuffle', true);

    switch(Method)
      case 3
        UbOutVname = '/u_basic';      
        VbOutVname = '/v_basic';      
  
        h5create(OutFile, UbOutVname, Vsize, 'DataType', 'single', 'ChunkSize', Csize, 'Deflate', 6, 'Shuffle', true);
        h5create(OutFile, VbOutVname, Vsize, 'DataType', 'single', 'ChunkSize', Csize, 'Deflate', 6, 'Shuffle', true);

        UdOutVname = '/u_dist';      
        VdOutVname = '/v_dist';      
  
        h5create(OutFile, UdOutVname, Vsize, 'DataType', 'single', 'ChunkSize', Csize, 'Deflate', 6, 'Shuffle', true);
        h5create(OutFile, VdOutVname, Vsize, 'DataType', 'single', 'ChunkSize', Csize, 'Deflate', 6, 'Shuffle', true);

        RadOutVname = '/radius';      
        h5create(OutFile, RadOutVname, Vsize, 'DataType', 'single', 'ChunkSize', Csize, 'Deflate', 6, 'Shuffle', true);

        AngOutVname = '/angle';      
        h5create(OutFile, AngOutVname, Vsize, 'DataType', 'single', 'ChunkSize', Csize, 'Deflate', 6, 'Shuffle', true);

        UavOutVname = '/u_vortex';      
        VavOutVname = '/v_vortex';      
        h5create(OutFile, UavOutVname, Vsize, 'DataType', 'single', 'ChunkSize', Csize, 'Deflate', 6, 'Shuffle', true);
        h5create(OutFile, VavOutVname, Vsize, 'DataType', 'single', 'ChunkSize', Csize, 'Deflate', 6, 'Shuffle', true);
    end

    fprintf('  Extracting vortex\n');
    for it = 1:Nt
      Start     = [ 1 1 1 it ];
      IoCount   = [ Nx Ny Nz 1 ];
      FiltCount = [ Nx Ny 1 1 ];

      U = squeeze(h5read(UinFile, Uvname, Start, IoCount));
      V = squeeze(h5read(VinFile, Vvname, Start, IoCount));

      FILTER = squeeze(h5read(FinFile, Fvname, Start, FiltCount));
      FILTER_3D = repmat(FILTER, [ 1 1 Nz ]);

      h5write(OutFile, UoutVname, U, Start, IoCount);
      h5write(OutFile, VoutVname, V, Start, IoCount);

      U_NV = single(zeros([ Nx Ny Nz ]));
      V_NV = single(zeros([ Nx Ny Nz ]));

      % Each method reads from U and V and places results in U_NV and V_NV.
      % Each method can write intermediate results as needed.
      switch(Method)
        case 1
          % Replace entries in U and V with nans wherever filter is a '1'
          U(FILTER_3D == 1) = nan;
          V(FILTER_3D == 1) = nan;

          % Interpolate level by level
          for k = 1:Nz
            U_NV(:,:,k) = single(inpaintn(squeeze(U(:,:,k))));
            V_NV(:,:,k) = single(inpaintn(squeeze(V(:,:,k))));
          end

        case 2
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          % WARNING - WORK IN PROGRESS
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          
          % Subtract the storm motion from U and V
          %U = U - StormSpeedX(it);
          %V = V - StormSpeedY(it);

          % Translate the cartesian grid origin to the storm center
          %    Reverse x and y arguments so that each row is a different x value and each
          %    column is a different y value. This keeps the size of X_GRID and Y_GRID consistent
          %    with the 3D data in the HDF5 file.
          StormX = X - StormLocX(it);
          StormY = Y - StormLocY(it);
          [ Y_GRID X_GRID ] = meshgrid(StormY, StormX);

          % Convert the cartesian grid to cylindrical
          %   THETA is now the angle from the positive x-axis (origin at storm center) to each grid cell
          [ THETA RADIUS ] = cart2pol(X_GRID, Y_GRID);
          THETA = repmat(THETA, [ 1 1 Nz ]);
          % note that THETA is (x,y,z) while RADIUS is (x,y)

          % PHI is the angle from the positive x axis to the wind vector
          PHI = atan2(V, U);

          % ALPHA is the angle to the wind vector relative to the radii from the storm center
          ALPHA = PHI - THETA;

          % Magnitude of wind vector
          WindMag = sqrt(U.^2 + V.^2);

          % tangential and radial winds
          VR = WindMag .* cos(ALPHA);
          VT = WindMag .* sin(ALPHA);

          % Azimuthal average
          Rstart = 0;
          Rinc = 1;
          Rend = 10;
          Rbins = Rstart:Rinc:Rend;
          Nb = length(Rbins);

          % Create a 2D array that shows which radial bands
          % each x,y location belongs to.
          R_SELECT = int32(RADIUS ./ Rinc) + 1;

          AzavgVr = zeros([ Nb Nz ]);
          AzavgVt = zeros([ Nb Nz ]);

%          for ib = 1:Nb-1
%            for iz = 1:Nz
%              for iy = 1:Ny
%                for ix = 1:Nx
%                end
%              end
%            end
%          end
%
%          for ir = 1:Nb-1
%            Select = (RADIUS >= Rbins(ir)) & (RADIUS < Rbins(ir+1));
%            for iz = 1:Nz
%              VR_SLICE = squeeze(VR(:,:,iz));
%              VT_SLICE = squeeze(VT(:,:,iz));
%              AzavgVr(ir,iz) = mean(VR_SLICE(Select));
%              AzavgVt(ir,iz) = mean(VT_SLICE(Select));
%            end
%          end

        case 3
          % Create the basic field by smoothing the input fields.
          [ U_B ] = KbrSmoothing(U,X,Y,Z);
          [ V_B ] = KbrSmoothing(V,X,Y,Z);
          h5write(OutFile, UbOutVname, U_B, Start, IoCount);
          h5write(OutFile, VbOutVname, V_B, Start, IoCount);

          % Disturbance field is original field minus the basic field
          U_D = U - U_B;
          V_D = V - V_B;
          h5write(OutFile, UdOutVname, U_D, Start, IoCount);
          h5write(OutFile, VdOutVname, V_D, Start, IoCount);

          % Filter out the hurricane piece from the disturbance field
          % to produce the non-hurricane piece of the disturbance field.

          % Radius from storm center
          StormX = X - StormLocX(it);
          StormY = Y - StormLocY(it);
          [ Y_GRID X_GRID ] = meshgrid(StormY, StormX);
          RADIUS = sqrt(X_GRID.^2 + Y_GRID.^2);
          RADIUS = repmat(RADIUS, [ 1 1 Nz ]);
          h5write(OutFile, RadOutVname, RADIUS, Start, IoCount);

          % Angle from storm center, round to nearest 10 degrees
          % Need to round angle to nearest 10 degrees to allow selection
          % of points along the perimeter of the storm (R0) in each sector
          ANGLE = atan2d(Y_GRID,X_GRID);
          A_SELECT = (ANGLE < 0);
          ANGLE(A_SELECT) = ANGLE(A_SELECT) + 360;
          ANGLE = round(ANGLE,-1); % round to nearest 10 degrees
          ANGLE = repmat(ANGLE, [ 1 1 Nz ]);
          h5write(OutFile, AngOutVname, ANGLE, Start, IoCount);

          % Eqn (3.8) in Kurihara et al. (1993)
          E = (exp(-(R0 - RADIUS).^2./L^2) - exp(-(R0^2)/L^2)) ./ (1 - exp(-(R0^2)/L^2));
          E(RADIUS > R0) = nan;

          % hD-bar quantities [Eqn (3.9) in Kurihara et al. (1993)]
          R_SELECT = (RADIUS >= R0*0.95) & (RADIUS <= R0*1.05);
          U_D_BAR = mean(U_D(R_SELECT));
          V_D_BAR = mean(V_D(R_SELECT));

          % hD(r0,theta)
          % Select along perimeter of storm in each 10 degree sector.
          % Take average of these selected points and fill in the entire
          % sector with that average value.
          U_D_R0 = single(zeros([ Nx Ny Nz ]));
          V_D_R0 = single(zeros([ Nx Ny Nz ]));
          for i = 0:10:360
            % peripheral values with angle equal to i degrees
            R_SELECT = (ANGLE == i) & (RADIUS >= R0*0.95) & (RADIUS <= R0*1.05);
            U_D_R0_VAL = mean(U_D(R_SELECT));
            V_D_R0_VAL = mean(V_D(R_SELECT));

            % Place the peripheral value in all locations where ANGLE == i
            A_SELECT = (ANGLE == i);
            U_D_R0(A_SELECT) = U_D_R0_VAL;
            V_D_R0(A_SELECT) = V_D_R0_VAL;
          end

          % Form analyzed vortex [Eqn (3.7) in Kurihara et al. (1993)]
          % E had nans outside the vortex periphery which will get
          % passed on to U_AV and V_AV.
          % Replace these with zeros so that U_AV and V_AV can be
          % subtracted from U and V.
          U_AV = U_D - (U_D_R0.*E + U_D_BAR.*(E-1));
          U_AV(isnan(U_AV)) = 0;
          V_AV = V_D - (V_D_R0.*E + V_D_BAR.*(E-1));
          V_AV(isnan(V_AV)) = 0;
          h5write(OutFile, UavOutVname, U_AV, Start, IoCount);
          h5write(OutFile, VavOutVname, V_AV, Start, IoCount);
          
          % Form environmental fields, U_NV and V_NV, by subtracting
          % the analyzed vortex from the original fields, U and V.
          U_NV = U - U_AV;
          V_NV = V - V_AV;

        otherwise
          if (it == 1)
            fprintf('WARNING: do not recognize method number %d, copying U,V to output\n', Method);
            fprintf('\n');
          end

          U_NV = U;
          V_NV = V;
      end

      h5write(OutFile, UnvOutVname, U_NV, Start, IoCount);
      h5write(OutFile, VnvOutVname, V_NV, Start, IoCount);

      if (mod(it, 10) == 0)
        fprintf('    Completed timestep: %d\n', it);
      end
    end
    fprintf('    Processed %d timesteps\n', it);
    fprintf('\n');

    % attach dimensions to output u,v for GRADS
    fprintf('  Writing %s\n', OutFile);

    DimOrder = { 'x' 'y' 'z' 't' };
    DimNames = 't z y x';

    fprintf('    %s\n', UoutVname);
    Units = 'm/s';
    LongName = 'Zonal Wind';
    AttachDimensionsXyzt(OutFile, UoutVname, DimOrder, Xvname, Yvname, Zvname, Tvname);
    NotateVariableXyzt(OutFile, UoutVname, Units, LongName, DimNames);

    fprintf('    %s\n', VoutVname);
    Units = 'm/s';
    LongName = 'Meridional Wind';
    AttachDimensionsXyzt(OutFile, VoutVname, DimOrder, Xvname, Yvname, Zvname, Tvname);
    NotateVariableXyzt(OutFile, VoutVname, Units, LongName, DimNames);

    fprintf('    %s\n', UnvOutVname);
    Units = 'm/s';
    LongName = 'Zonal Wind No Vortex';
    AttachDimensionsXyzt(OutFile, UnvOutVname, DimOrder, Xvname, Yvname, Zvname, Tvname);
    NotateVariableXyzt(OutFile, UnvOutVname, Units, LongName, DimNames);

    fprintf('    %s\n', VnvOutVname);
    Units = 'm/s';
    LongName = 'Meridional Wind No Vortex';
    AttachDimensionsXyzt(OutFile, VnvOutVname, DimOrder, Xvname, Yvname, Zvname, Tvname);
    NotateVariableXyzt(OutFile, VnvOutVname, Units, LongName, DimNames);

    switch(Method)
      case 3
        fprintf('    %s\n', UbOutVname);
        Units = 'm/s';
        LongName = 'Zonal Wind Smoothed';
        AttachDimensionsXyzt(OutFile, UbOutVname, DimOrder, Xvname, Yvname, Zvname, Tvname);
        NotateVariableXyzt(OutFile, UbOutVname, Units, LongName, DimNames);

        fprintf('    %s\n', VbOutVname);
        Units = 'm/s';
        LongName = 'Meridional Wind Smoothed';
        AttachDimensionsXyzt(OutFile, VbOutVname, DimOrder, Xvname, Yvname, Zvname, Tvname);
        NotateVariableXyzt(OutFile, VbOutVname, Units, LongName, DimNames);

        fprintf('    %s\n', UdOutVname);
        Units = 'm/s';
        LongName = 'Zonal Wind Perturbation';
        AttachDimensionsXyzt(OutFile, UdOutVname, DimOrder, Xvname, Yvname, Zvname, Tvname);
        NotateVariableXyzt(OutFile, UdOutVname, Units, LongName, DimNames);

        fprintf('    %s\n', VdOutVname);
        Units = 'm/s';
        LongName = 'Meridional Wind Perturbation';
        AttachDimensionsXyzt(OutFile, VdOutVname, DimOrder, Xvname, Yvname, Zvname, Tvname);
        NotateVariableXyzt(OutFile, VdOutVname, Units, LongName, DimNames);

        fprintf('    %s\n', RadOutVname);
        Units = 'm';
        LongName = 'Radius relative to storm center';
        AttachDimensionsXyzt(OutFile, RadOutVname, DimOrder, Xvname, Yvname, Zvname, Tvname);
        NotateVariableXyzt(OutFile, RadOutVname, Units, LongName, DimNames);

        fprintf('    %s\n', AngOutVname);
        Units = 'degrees';
        LongName = 'Angle relative to storm center';
        AttachDimensionsXyzt(OutFile, AngOutVname, DimOrder, Xvname, Yvname, Zvname, Tvname);
        NotateVariableXyzt(OutFile, AngOutVname, Units, LongName, DimNames);

        fprintf('    %s\n', UavOutVname);
        Units = 'm/s';
        LongName = 'Zonal Wind No Vortex';
        AttachDimensionsXyzt(OutFile, UavOutVname, DimOrder, Xvname, Yvname, Zvname, Tvname);
        NotateVariableXyzt(OutFile, UavOutVname, Units, LongName, DimNames);

        fprintf('    %s\n', VavOutVname);
        Units = 'm/s';
        LongName = 'Meridional Wind No Vortex';
        AttachDimensionsXyzt(OutFile, VavOutVname, DimOrder, Xvname, Yvname, Zvname, Tvname);
        NotateVariableXyzt(OutFile, VavOutVname, Units, LongName, DimNames);
    end

  end % cases
end % function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% KbrSmoothing()
%
% This function performs the smoothing described in Kurihara et al. (1993)
% on a single field.
%
function [ VarSmooth ] = KbrSmoothing(Var,X,Y,Z)
  % K coefficient values from Kurihara et al. (1993)
  M = [ 2 3 4 2 5 6 7 2 8 9 2 ];
  K = 0.5 ./ (1-cos((2*pi)./M));
  Nk = length(K);

  % Kurihara et al. (1993) specified the formula for 1 degree grid spacing, and for
  % filtering out components up to 9 degrees wavelength. This data is on a 1/3 degree
  % grid spacing. When I tried to adjust the formula so that it would still filter
  % up to 9 degree wavelength, the function blew up. To deal with this, sample the
  % input Var every 3rd point to get gridded data on a 1 degree grid.
  Vinc = 3;
  VarSmoothLow = Var(1:Vinc:end,1:Vinc:end,:);

  % First smooth in the x direction
  for i = 1:Nk
    VarIm1 = VarSmoothLow(1:end-2,:,:);
    VarI   = VarSmoothLow(2:end-1,:,:);
    VarIp1 = VarSmoothLow(3:end,  :,:);

    VarSmoothLow(2:end-1,:,:) = VarSmoothLow(2:end-1,:,:) + K(i) .* (VarIm1 + VarIp1 - 2.*VarI);
  end
  clear VarIm1;
  clear VarI;
  clear VarIp1;

  % Second smooth in the y direction
  for i = 1:Nk
    VarJm1 = VarSmoothLow(:,1:end-2,:);
    VarJ   = VarSmoothLow(:,2:end-1,:);
    VarJp1 = VarSmoothLow(:,3:end,  :);

    VarSmoothLow(:,2:end-1,:) = VarSmoothLow(:,2:end-1,:) + K(i) .* (VarJm1 + VarJp1 - 2.*VarJ);
  end
  clear VarJm1;
  clear VarJ;
  clear VarJp1;

  % VarSmoothLow is at 1/3 resolution than Var. Interpolate back to the higher
  % resolution so that it matches Var. [XYZ]L are gridded coordinate values
  % for low resolution grid, [XYZ]H are same for high resolution grid.
  [ XL, YL, ZL ] = ndgrid(X(1:Vinc:end),Y(1:Vinc:end),Z);
  [ XH, YH, ZH ] = ndgrid(X,Y,Z);
  VarSmooth = interpn(XL,YL,ZL,VarSmoothLow,XH,YH,ZH);
end
