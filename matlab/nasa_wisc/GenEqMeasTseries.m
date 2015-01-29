function [ ] = GenEqMeasTseries(ConfigFile)
% GenEqMeasTseries generate time series of equllibrium measurements

  % Read the config file to get the structure of how the data is laid out in
  % the file system.
  [ Config ] = ReadConfig(ConfigFile);

  Ddir = Config.DiagDir;
  Hdir = 'HDF5';

  % make sure output directory exists
  if (exist(Ddir, 'dir') ~= 7)
    mkdir(Ddir);
  end

  % cases
  %CaseList = { 'RCE50_RECT' };
  CaseList = { 'RCE50_RECT_S303' };
  %CaseList = { 'RCE50_SQ' };
  Nc = length(CaseList);

  % input file specs
  % thermal heat flux
  SfcLatFbase  = 'lat_flux-a-AS-2012-01-01-000000-g1.h5';
  SfcLatVname  = 'lat_flux';
  SfcSensFbase = 'sens_flux-a-AS-2012-01-01-000000-g1.h5';
  SfcSensVname = 'sens_flux';

  % radiative flux divergence
  SfcSwdnFbase = 'rshort-a-AS-2012-01-01-000000-g1.h5';
  SfcSwdnVname = 'rshort';
  SfcLwdnFbase = 'rlong-a-AS-2012-01-01-000000-g1.h5';
  SfcLwdnVname = 'rlong';
  SfcLwupFbase = 'rlongup-a-AS-2012-01-01-000000-g1.h5';
  SfcLwupVname = 'rlongup';
  SfcAlbFbase  = 'albedt-a-AS-2012-01-01-000000-g1.h5';
  SfcAlbVname  = 'albedt';

  TopSwdnFbase = 'swdn_toa-a-AS-2012-01-01-000000-g1.h5';
  TopSwdnVname = 'swdn';
  TopSwupFbase = 'swup_toa-a-AS-2012-01-01-000000-g1.h5';
  TopSwupVname = 'swup';
  TopLwupFbase = 'lwup_toa-a-AS-2012-01-01-000000-g1.h5';
  TopLwupVname = 'lwup';

  % output file specs
  OutFprefix = 'eq_meas';
  ThfVname = 'therm_heat_flux';
  QradVname = 'rad_flux_div';

  fprintf('***************************************************************\n');
  fprintf('Generating domain average time series:\n');

  for icase = 1:Nc
    Case = CaseList{icase};
    OutFile = sprintf('%s/%s_%s.h5', Ddir, OutFprefix, Case);

    fprintf('  Case: %s\n', Case);
    fprintf('\n');

    % Assemble input file names
    SfcLatFile  = sprintf('%s/%s/%s', Hdir, Case, SfcLatFbase);
    SfcSensFile = sprintf('%s/%s/%s', Hdir, Case, SfcSensFbase);

    SfcSwdnFile = sprintf('%s/%s/%s', Hdir, Case, SfcSwdnFbase);
    SfcLwdnFile = sprintf('%s/%s/%s', Hdir, Case, SfcLwdnFbase);
    SfcLwupFile = sprintf('%s/%s/%s', Hdir, Case, SfcLwupFbase);
    SfcAlbFile  = sprintf('%s/%s/%s', Hdir, Case, SfcAlbFbase);

    TopSwdnFile  = sprintf('%s/%s/%s', Hdir, Case, TopSwdnFbase);
    TopSwupFile  = sprintf('%s/%s/%s', Hdir, Case, TopSwupFbase);
    TopLwupFile  = sprintf('%s/%s/%s', Hdir, Case, TopLwupFbase);

    fprintf('  Reading input files:\n');
    fprintf('    %s (%s)\n', SfcLatFile, SfcLatVname);
    fprintf('    %s (%s)\n', SfcSensFile, SfcSensVname);

    fprintf('    %s (%s)\n', SfcSwdnFile, SfcSwdnVname);
    fprintf('    %s (%s)\n', SfcLwdnFile, SfcLwdnVname);
    fprintf('    %s (%s)\n', SfcLwupFile, SfcLwupVname);
    fprintf('    %s (%s)\n', SfcAlbFile, SfcAlbVname);

    fprintf('    %s (%s)\n', TopSwdnFile, TopSwdnVname);
    fprintf('    %s (%s)\n', TopSwupFile, TopSwupVname);
    fprintf('    %s (%s)\n', TopLwupFile, TopLwupVname);
    fprintf('\n');

    % Set up to read input files one time step at a time
    SFC_LAT_DS  = ncgeodataset(SfcLatFile);
    SFC_SENS_DS = ncgeodataset(SfcSensFile);

    SFC_SWDN_DS = ncgeodataset(SfcSwdnFile);
    SFC_LWDN_DS = ncgeodataset(SfcLwdnFile);
    SFC_LWUP_DS = ncgeodataset(SfcLwupFile);
    SFC_ALB_DS  = ncgeodataset(SfcAlbFile);
    TOP_SWDN_DS = ncgeodataset(TopSwdnFile);
    TOP_SWUP_DS = ncgeodataset(TopSwupFile);
    TOP_LWUP_DS = ncgeodataset(TopLwupFile);


    SFC_LAT_VAR  = SFC_LAT_DS.geovariable(SfcLatVname);
    SFC_SENS_VAR = SFC_SENS_DS.geovariable(SfcSensVname);

    SFC_SWDN_VAR = SFC_SWDN_DS.geovariable(SfcSwdnVname);
    SFC_LWDN_VAR = SFC_LWDN_DS.geovariable(SfcLwdnVname);
    SFC_LWUP_VAR = SFC_LWUP_DS.geovariable(SfcLwupVname);
    SFC_ALB_VAR  = SFC_ALB_DS.geovariable(SfcAlbVname);
    TOP_SWDN_VAR = TOP_SWDN_DS.geovariable(TopSwdnVname);
    TOP_SWUP_VAR = TOP_SWUP_DS.geovariable(TopSwupVname);
    TOP_LWUP_VAR = TOP_LWUP_DS.geovariable(TopLwupVname);

    T_VAR = SFC_LAT_DS.geovariable('/t_coords');


    % First read in the time coordinates and determine the number of time steps
    T = T_VAR.data(:) ./ (3600 * 24);  % convert sec to days
    Nt = length(T);

    % For each time step calculate the domain average of 
    % Variables are either (t,z,y,x) or (t,y,x) where z is size 1
    THF  = zeros([ Nt 1 ]);
    QRAD = zeros([ Nt 1 ]);

    for it = 1:Nt
      % Read in the variables -> after squeeze will be (y,x)
      SFC_LAT  = squeeze(SFC_LAT_VAR.data(it,:,:));
      SFC_SENS = squeeze(SFC_SENS_VAR.data(it,:,:));

      SFC_SWDN = squeeze(SFC_SWDN_VAR.data(it,:,:));
      SFC_LWDN = squeeze(SFC_LWDN_VAR.data(it,:,:));
      SFC_LWUP = squeeze(SFC_LWUP_VAR.data(it,:,:));
      SFC_ALB  = squeeze(SFC_ALB_VAR.data(it,:,:));
      TOP_SWDN = squeeze(TOP_SWDN_VAR.data(it,:,:,:));
      TOP_SWUP = squeeze(TOP_SWUP_VAR.data(it,:,:,:));
      TOP_LWUP = squeeze(TOP_LWUP_VAR.data(it,:,:,:));

      % Form the THF and QRAD numbers
      %  THF is the sum of the surface lat heat flux and sensible heat flux
      %
      %  QRAD is the net radiative flux of the entire column
      %    or: (Net upward flux at TOA) - (Net upward flux at surface)
      %
      %    where:
      %      (Net upward flux at TOA) = (TOP_SWUP + TOP_LWUP) - (TOP_SWDN)
      %      (Net upward flux at surface) = (SFC_SWUP + SFC_LWUP) - (SFC_SWDN + SFC_LWDN)
      %
      %    where:
      %      SFC_SWUP = SFC_SWDN * SFC_ALB
      %
      % Calculate THF, QRAD per column and then take average
      %
      TempVar = SFC_LAT + SFC_SENS;
      THF(it) = mean(TempVar(:));

      TempVar = ((TOP_SWUP + TOP_LWUP) - (TOP_SWDN)) - (((SFC_SWDN.*SFC_ALB) + SFC_LWUP) - (SFC_SWDN + SFC_LWDN));
      QRAD(it) = mean(TempVar(:));

      if (mod(it,50) == 0)
        fprintf ('  Completed timestep %d out of %d\n', it, Nt);
      end
    end
    fprintf('\n');

    % Write out dummy coordinate values to keep ReadSelectXyzt happy
    Xdummy = 1;
    Ydummy = 1;
    Zdummy = 1;

    fprintf('    Writing: %s (%s, %s)\n', OutFile, ThfVname, QradVname);
    fprintf('\n');

    OutVar = reshape(THF, [ 1 1 1 Nt ]);
    hdf5write(OutFile, ThfVname, OutVar); 

    OutVar = reshape(QRAD, [ 1 1 1 Nt ]);
    hdf5write(OutFile, QradVname, OutVar, 'WriteMode', 'append'); 

    hdf5write(OutFile, '/x_coords', Xdummy, 'WriteMode', 'append');
    hdf5write(OutFile, '/y_coords', Ydummy, 'WriteMode', 'append');
    hdf5write(OutFile, '/z_coords', Zdummy, 'WriteMode', 'append');
    hdf5write(OutFile, '/t_coords', T,      'WriteMode', 'append');
  end
end