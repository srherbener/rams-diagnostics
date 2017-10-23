function [ ] = GenPtrackMomentum()

  CaseList = {
   'TSD_SAL_DUST'
   'TSD_SAL_NODUST'
   'TSD_NONSAL_DUST'
   'TSD_NONSAL_NODUST'
   };
  
  Nc = length(CaseList);

  PreSalStart = 10;  % start, end Pre-SAL time period, sim time in hours
  PreSalEnd   = 30;
  SalStart    = 40;  % start, end SAL time period
  SalEnd      = 60;

  Xname = '/x_coords';
  Yname = '/y_coords';
  Zname = '/z_coords';
  Tname = '/t_coords';

  Pdx = 1;
  Pdy = 2.593;
  Pslope = Pdy / Pdx;
  Pangle = atan2(Pdy, Pdx);
  
  fprintf('Translating horizontal velocity vectors to ptrack axes\n');
  fprintf('  Ptrack dx, dy: %f %f\n', Pdx, Pdy);
  fprintf('  Ptrack slope: %f\n', Pslope);
  fprintf('  Ptrack angle: %f\n', Pangle);
  fprintf('\n');

  % Loop through all time steps
  %   Read in u, v (which are on the domain axes)
  %   Convert to polar form relative to domain axes
  %   Translate angle to be relative to ptrack axes
  %   Convert back to cartesian form (which will now be on ptrack axes)
  %
  for icase = 1:Nc
    Case = CaseList{icase};
    fprintf('***************************************\n');
    fprintf('Case: %s\n', Case);

    %Ufname = sprintf('XsectionData/ptrack_u_%s.h5', Case);
    Ufname = sprintf('XsectionData/ptrack_u_p_%s.h5', Case);
    Uvname = '/u';

    %Vfname = sprintf('XsectionData/ptrack_v_%s.h5', Case);
    Vfname = sprintf('XsectionData/ptrack_v_p_%s.h5', Case);
    Vvname = '/v';

    fprintf('  Reading: %s (%s)\n', Ufname, Uvname);
    fprintf('  Reading: %s (%s)\n', Vfname, Vvname);
    fprintf('\n');

    % Read in coordinates from u file (v will have same coords)
    X = squeeze(h5read(Ufname, Xname));
    Y = squeeze(h5read(Ufname, Yname));
    Z = squeeze(h5read(Ufname, Zname));
    T = squeeze(h5read(Ufname, Tname));

    Nx = length(X);
    Ny = length(Y);
    Nz = length(Z);
    Nt = length(T);

    % Find start,end indices for time periods
    SIM_T = T ./ 3600 - 42;
    PS_T1 = find(SIM_T >= PreSalStart, 1, 'first');
    PS_T2 = find(SIM_T <= PreSalEnd,   1, 'last');

    S_T1 = find(SIM_T >= SalStart, 1, 'first');
    S_T2 = find(SIM_T <= SalEnd,   1, 'last');

    % Read in and translate from domain axes to ptrack axes
    U_DOM = squeeze(h5read(Ufname, Uvname));
    V_DOM = squeeze(h5read(Vfname, Vvname));

    % Convert to polar coords (relative to domain axes)
    MAG   = sqrt(U_DOM.^2 + V_DOM.^2);
    ANGLE = atan2(V_DOM, U_DOM);

    % Translate to ptrack coords
    ANGLE = ANGLE - Pangle;

    % Convert to cartesian coords (relative to ptrack axes)
    U_PTRACK = MAG .* cos(ANGLE);
    V_PTRACK = MAG .* sin(ANGLE);

    % Create temporal averages
    PS_U_PTRACK = squeeze(mean(U_PTRACK(:,:,PS_T1:PS_T2),3));
    PS_V_PTRACK = squeeze(mean(V_PTRACK(:,:,PS_T1:PS_T2),3));

    S_U_PTRACK = squeeze(mean(U_PTRACK(:,:,S_T1:S_T2),3));
    S_V_PTRACK = squeeze(mean(V_PTRACK(:,:,S_T1:S_T2),3));

    % Output, get rid of the dummy y dimension
    %OutFname = sprintf('DIAGS/ptrack_hvelocity_%s.h5', Case);
    OutFname = sprintf('DIAGS/ptrack_hvelocity_p_%s.h5', Case);

    % Remove output file if it exists so that new dataset can be written.
    if (exist(OutFname, 'file') == 2)
      delete(OutFname);
    end

    % Write coordinates into output file
    CreateDimensionsXyzt(OutFname, X, Y, Z, T, Xname, Yname, Zname, Tname);
    NotateDimensionsXyzt(OutFname, Xname, Yname, Zname, Tname);

    % Time series
    Vsize = [ Nx Nz Nt ];
    DimOrder = { 'x' 'z' 't' };

    OutVname = '/u';
    fprintf('  Writing: %s (%s)\n', OutFname, OutVname);
    h5create(OutFname, OutVname, Vsize)
    h5write(OutFname, OutVname, U_PTRACK);
    AttachDimensionsXyzt(OutFname, OutVname, DimOrder, Xname, Yname, Zname, Tname);

    OutVname = '/v';
    fprintf('  Writing: %s (%s)\n', OutFname, OutVname);
    h5create(OutFname, OutVname, Vsize)
    h5write(OutFname, OutVname, V_PTRACK);
    AttachDimensionsXyzt(OutFname, OutVname, DimOrder, Xname, Yname, Zname, Tname);

    % Temporal Averages
    Vsize = [ Nx Nz ];
    DimOrder = { 'x' 'z' };

    OutVname = '/ps_u';
    fprintf('  Writing: %s (%s)\n', OutFname, OutVname);
    h5create(OutFname, OutVname, Vsize)
    h5write(OutFname, OutVname, PS_U_PTRACK);
    AttachDimensionsXyzt(OutFname, OutVname, DimOrder, Xname, Yname, Zname, Tname);

    OutVname = '/ps_v';
    fprintf('  Writing: %s (%s)\n', OutFname, OutVname);
    h5create(OutFname, OutVname, Vsize)
    h5write(OutFname, OutVname, PS_V_PTRACK);
    AttachDimensionsXyzt(OutFname, OutVname, DimOrder, Xname, Yname, Zname, Tname);

    OutVname = '/s_u';
    fprintf('  Writing: %s (%s)\n', OutFname, OutVname);
    h5create(OutFname, OutVname, Vsize)
    h5write(OutFname, OutVname, S_U_PTRACK);
    AttachDimensionsXyzt(OutFname, OutVname, DimOrder, Xname, Yname, Zname, Tname);

    OutVname = '/s_v';
    fprintf('  Writing: %s (%s)\n', OutFname, OutVname);
    h5create(OutFname, OutVname, Vsize)
    h5write(OutFname, OutVname, S_V_PTRACK);
    AttachDimensionsXyzt(OutFname, OutVname, DimOrder, Xname, Yname, Zname, Tname);

    fprintf('\n');
  end
end
