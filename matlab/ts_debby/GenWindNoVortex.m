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
%         a) Kurihara et al. (1995) says to use sigma = 0.85 level
%            which is rougly 1500 m elevation. Say p0 = 1000 mb, then
%            sigma = 0.85 is 850 mb which is roughly 1500 m elevation

  Method = 3;

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
    UnvOutVname = '/u_nv';
    VnvOutVname = '/v_nv';

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

    % Read in the storm speed components
    if (Method == 2)
      StormSpeedX = squeeze(h5read(ScFile, SsxVname));
      StormSpeedY = squeeze(h5read(ScFile, SsyVname));

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
    fprintf('    %s\n', UoutVname);
    fprintf('    %s\n', VoutVname);
    fprintf('    %s\n', UnvOutVname);
    fprintf('    %s\n', VnvOutVname);
    AttachDimensionsXyzt(OutFile, UoutVname, DimOrder, Xvname, Yvname, Zvname, Tvname);
    AttachDimensionsXyzt(OutFile, VoutVname, DimOrder, Xvname, Yvname, Zvname, Tvname);
    AttachDimensionsXyzt(OutFile, UnvOutVname, DimOrder, Xvname, Yvname, Zvname, Tvname);
    AttachDimensionsXyzt(OutFile, VnvOutVname, DimOrder, Xvname, Yvname, Zvname, Tvname);

    switch(Method)
      case 3
        fprintf('    %s\n', UbOutVname);
        fprintf('    %s\n', VbOutVname);
        AttachDimensionsXyzt(OutFile, UbOutVname, DimOrder, Xvname, Yvname, Zvname, Tvname);
        AttachDimensionsXyzt(OutFile, VbOutVname, DimOrder, Xvname, Yvname, Zvname, Tvname);

        fprintf('    %s\n', UdOutVname);
        fprintf('    %s\n', VdOutVname);
        AttachDimensionsXyzt(OutFile, UdOutVname, DimOrder, Xvname, Yvname, Zvname, Tvname);
        AttachDimensionsXyzt(OutFile, VdOutVname, DimOrder, Xvname, Yvname, Zvname, Tvname);
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
