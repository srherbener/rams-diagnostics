function [ ] = GenPtrackMomentum()

  CaseList = {
   'TSD_SAL_DUST'
   'TSD_SAL_NODUST'
   'TSD_NONSAL_DUST'
   'TSD_NONSAL_NODUST'
   };
  
  Nc = length(CaseList);

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

    Ufname = sprintf('XsectionData/ptrack_u_%s.h5', Case);
    Uvname = '/u';

    Vfname = sprintf('XsectionData/ptrack_v_%s.h5', Case);
    Vvname = '/v';

    OutFname = sprintf('DIAGS/ptrack_hvelocity_%s.h5', Case);
    OutUvname = '/u';
    OutVvname = '/v';

    fprintf('  Reading: %s (%s)\n', Ufname, Uvname);
    fprintf('  Reading: %s (%s)\n', Vfname, Vvname);
    fprintf('  Writing: %s (%s, %s)\n', OutFname, OutUvname, OutVvname);
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

    % Remove output file if it exists so that new dataset can be written.
    if (exist(OutFname, 'file') == 2)
      delete(OutFname);
    end

    % Write coordinates into output file
    CreateDimensionsXyzt(OutFname, X, Y, Z, T, Xname, Yname, Zname, Tname);
    NotateDimensionsXyzt(OutFname, Xname, Yname, Zname, Tname);

    % Set up datasets that can have a 3D field written one time step at a time.
    Vsize = [ Nx Nz Inf ];
    Csize = [ Nx Nz 1 ];
    h5create(OutFname, OutUvname, Vsize, 'ChunkSize', Csize, 'Deflate', 6, 'Shuffle', true);
    h5create(OutFname, OutVvname, Vsize, 'ChunkSize', Csize, 'Deflate', 6, 'Shuffle', true);

    % For output, y is dummy coordinate so remove it in the output file
    InCount  = [ Nx Ny Nz 1 ];
    OutCount = [ Nx Nz 1 ];
    
    % Read in x-z plane, translate to ptrack coords
    for it = 1:Nt
      InStart = [ 1 1 1 it ];
      U_DOM = squeeze(h5read(Ufname, Uvname, InStart, InCount));
      V_DOM = squeeze(h5read(Vfname, Vvname, InStart, InCount));

      % Convert to polar coords (relative to domain axes)
      MAG   = sqrt(U_DOM.^2 + V_DOM.^2);
      ANGLE = atan2(V_DOM, U_DOM);

      % Translate to ptrack coords
      ANGLE = ANGLE - Pangle;

      % Convert to cartesian coords (relative to ptrack axes)
      U_PTRACK = MAG .* cos(ANGLE);
      V_PTRACK = MAG .* sin(ANGLE);

      % Output is x-z plane
      OutStart = [ 1 1 it ];
      h5write(OutFname, OutUvname, U_PTRACK, OutStart, OutCount);
      h5write(OutFname, OutVvname, V_PTRACK, OutStart, OutCount);

      if (mod(it,10) == 0)
        fprintf('  Completed time step: %d\n', it);
      end
    end
    fprintf('  Total time steps processed: %d\n', it);

    % Attach dimensions (COARDS format)
    DimOrder = { 'x' 'z' 't' };
    AttachDimensionsXyzt(OutFname, OutUvname, DimOrder, Xname, Yname, Zname, Tname);
    AttachDimensionsXyzt(OutFname, OutVvname, DimOrder, Xname, Yname, Zname, Tname);

    fprintf('\n');
  end
end
