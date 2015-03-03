function [ ] = GenNonSalRamsInfo()
% GenNonSalRamsInfo function to generate RAMS TS Debby non-SAL simulation data (prep for GenNonSalRamsInit.m)

  % Read in RAMS var files (dprep data analyzed to RAMS grids) and convert
  % the RV values to RH. Then mask off the region that was used in
  % GenNonSalDpFiles() to create the input for the RAMS MAKEVFILE run. Output
  % only the RH and THETA values, with RAMS (i,j,k) location, for RAMS to read in and
  % use for initialization.

  Ngrids = 3;
  Hdir = 'HDF5';
  Ddir = 'DP_FILES/NONSAL';

  OutFile = sprintf('%s/NonSalRamsInfo.h5', Hdir);
  
  fprintf('*****************************************************************************\n');
  fprintf('Extracting RH and THETA for RAMS TS Debby non-SAL simulation initialization:\n');
  fprintf('\n');


  % Read in the orignal grid and set up to do interpolation to the RAMS grids
  MaskFile = sprintf('%s/dp-p2006-08-22-0600-MASK.h5', Ddir);
  MaskVar = 'MASK';
  MaskLonVar = 'LON';
  MaskLatVar = 'LAT';
  MaskPressVar = 'PRESS';

  fprintf('  Reading dprep mask file: %s -> %s, %s, %s, %s\n', MaskFile, MaskVar, MaskLonVar, MaskLatVar, MaskPressVar);
  fprintf('\n');

  MASK = squeeze(hdf5read(MaskFile, MaskVar));
  MASK_LON = squeeze(hdf5read(MaskFile, MaskLonVar));
  MASK_LAT = squeeze(hdf5read(MaskFile, MaskLatVar));
  MASK_PRESS = squeeze(hdf5read(MaskFile, MaskPressVar));

  % The MASK grid is (Lon, Lat, Pressure) and all other grids are (Lon, Lat, Height). Need to
  % convert the pressure levels in the mask grid to height levels.
  MASK_Z = CalcAltitude(MASK_PRESS);

  % Create grid meshes for subsequent call to interpn
  [ MX, MY, MZ ] = ndgrid(MASK_LON, MASK_LAT, MASK_Z);

  % Write the original mask out for comparison purposes. If you get rid of this section, make sure
  % you write something into the output file so that the appending in the following loop doesn't
  % just keep adding more a more data to an existing file. Writing a header string while not using
  % the append write mode will be sufficient (ie, remove the existing file and create a new file).
  fprintf('  Writing orignal mask to: %s\n', OutFile);
  fprintf('\n');

  hdf5write(OutFile, 'ORIG_MASK',  MASK);
  hdf5write(OutFile, 'ORIG_LON',   MASK_LON,   'WriteMode', 'append');
  hdf5write(OutFile, 'ORIG_LAT',   MASK_LAT,   'WriteMode', 'append');
  hdf5write(OutFile, 'ORIG_PRESS', MASK_PRESS, 'WriteMode', 'append');
  hdf5write(OutFile, 'ORIG_Z',     MASK_Z,     'WriteMode', 'append');

  % Read in each grid and interpolate the mask data to that grid. Also, calculate the new RH
  % values from the data in the var files. Write out the interpolated masks and RH values for
  % processing into initialization data for RAMS.
  for igrid = 1:Ngrids
    GridFile = sprintf('%s/grid%d-TSD_DRY_NODUST-AS-2006-08-20-120000-g%d.h5', Hdir, igrid, igrid);
    RamsLonVar = 'x_coords';
    RamsLatVar = 'y_coords';
    RamsZvar = 'z_coords';

    VarFile = sprintf('%s/var-V-2006-08-22-060000-g%d-NONSAL.h5', Hdir, igrid);
    PiVar = 'PI';
    RvVar = 'RV';
    ThetaVar = 'THETA';

    fprintf('  Grid: %d\n', igrid);
    fprintf('    Reading grid file: %s -> %s, %s, %s\n', GridFile, RamsLonVar, RamsLatVar, RamsZvar);
    fprintf('    Reading var file: %s -> %s, %s, %s\n', VarFile, PiVar, RvVar, ThetaVar);
    fprintf('\n');

    RAMS_LON = squeeze(hdf5read(GridFile, RamsLonVar));
    RAMS_LAT = squeeze(hdf5read(GridFile, RamsLatVar));
    RAMS_Z = squeeze(hdf5read(GridFile, RamsZvar));
    Nx = length(RAMS_LON);
    Ny = length(RAMS_LAT);
    Nz = length(RAMS_Z);

    % The VAR files have the dimension descriptions reversed (RAMS "feature"), however
    % the data is in the correct order. Fix the dimensions by doing a reshape, in Nx, Ny, Nz order.
    RAMS_PI = squeeze(hdf5read(VarFile, PiVar));
    RAMS_RV = squeeze(hdf5read(VarFile, RvVar));
    RAMS_THETA = squeeze(hdf5read(VarFile, ThetaVar));
    RAMS_PI = reshape(RAMS_PI, [ Nx Ny Nz ]);
    RAMS_RV = reshape(RAMS_RV, [ Nx Ny Nz ]);
    RAMS_THETA = reshape(RAMS_THETA, [ Nx Ny Nz ]);

    % Interpolate the MASK data to the RAMS grid. Use interp3d.
    [ RX, RY, RZ ] = ndgrid(RAMS_LON, RAMS_LAT, RAMS_Z);
    RAMS_MASK = interpn(MX, MY, MZ, MASK, RX, RY, RZ);

    % The interpolation will cause regions to have values between 0 and 1. Turn this
    % back into zeros and ones using a threshold. This will also tend to push back
    % regions around the perimiter of the original boundary between zeros and ones
    % that spread outward (from the ones) due to the interpolation.
    RAMS_MASK = single(RAMS_MASK >= 0.75);

    % Convert the RV data into relative humidity
    RAMS_RH = CalcRhumidity(RAMS_RV, RAMS_PI, RAMS_THETA);

    % Output
    fprintf('  Writing grid %d to: %s\n', igrid, OutFile);
    fprintf('\n');
    
    OutVar = sprintf('%s_G%d', 'RAMS_MASK', igrid);
    hdf5write(OutFile, OutVar, RAMS_MASK, 'WriteMode', 'append');
    OutVar = sprintf('%s_G%d', 'RAMS_LON', igrid);
    hdf5write(OutFile, OutVar, RAMS_LON, 'WriteMode', 'append');
    OutVar = sprintf('%s_G%d', 'RAMS_LAT', igrid);
    hdf5write(OutFile, OutVar, RAMS_LAT, 'WriteMode', 'append');
    OutVar = sprintf('%s_G%d', 'RAMS_Z', igrid);
    hdf5write(OutFile, OutVar, RAMS_Z, 'WriteMode', 'append');

    OutVar = sprintf('%s_G%d', 'RAMS_RH', igrid);
    hdf5write(OutFile, OutVar, RAMS_RH, 'WriteMode', 'append');

    OutVar = sprintf('%s_G%d', 'RAMS_THETA', igrid);
    hdf5write(OutFile, OutVar, RAMS_THETA, 'WriteMode', 'append');
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Altimitry Eqauation to convert pressure to height. Based on
% hydrostatic balance and a given uniform lapse rate.
%
function [ Z ] = CalcAltitude(P)
% CalcAltitude - generate height from a given pressure, uses altimitry equation

  % Altimitry Eqn:
  %
  %   Z = (T0 / LapseRate) * [ 1 - (p/p0)^(R*LapseRate/g) ]
  %
  %      g = 9.81 m/s^2
  %     T0 = 288 K  (temp at sea level)
  %     LapseRate = 6.5 K/km   (~ psuedoadiabitic rate)
  %     R = dry air gas constant 287 J/kg/km
  %     p0 = 1013 mb  (press at sea level)
  %
  %   R * LapseRate / g = 0.1902
  %   T0 / LapseRate = 44310 m
  %
  % Formula becomes:
  %
  %   Z (in m) = 44310 * [ 1 - (p/1013)^(0.1902)]

  Z = 44310 .* (1 - (P./1013).^(0.1902));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% conversion from RV to RH
%
function [ RH ] = CalcRhumidity(RV, RAMS_PI, THETA)
% CalcRhumidity - generate relative humidity from given vapor mixing ratio, pi and theta

  % Constants
  Cp = 1004;               % heat capacity of dry air, J/kg/K
  Rd = 287;                % dry air gas constant, J/kg/K
  p0 = 1013;               % reference pressure, mb

  OneOverCp = 1 / Cp;
  CpOverRd = Cp / Rd;

  % RAMS uses PI*Cp
  PI = RAMS_PI .* OneOverCp;

  % Get P and T from PI and THETA
  %
  %  P = p0 * PI^(Cp/Rd)
  %  T = THETA * PI
  %
  % The 1/Cp factor is because RAMS uses PI*Cp
  P = p0 .* PI.^(CpOverRd);    % P is in mb
  T = THETA .* PI;             % T is in Kelvin
  T = T - 273.15;              % T is in Celsius

  % Saturation vapor pressure for given T (in Celsius)
  ES =  6.112 .* exp(17.67 .* T ./ (T + 243.5)); % ES is in mb

  % Saturation vapor mixing ratio (pressures in mb)
  RVS = 0.622 .* (ES ./ (P - ES)); % in kg/kg

  % Relative Humidity
  RH = RV ./ RVS; % leave in fractional form 
end
