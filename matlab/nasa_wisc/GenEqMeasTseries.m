function [ ] = GenEqMeasTseries()
% GenEqMeasTseries generate time series of equllibrium measurements

  Ddir = 'DIAGS';
  Hdir = 'HDF5';

  % make sure output directory exists
  if (exist(Ddir, 'dir') ~= 7)
    mkdir(Ddir);
  end

  % CaseList names each case and indicates where to find the
  % input variable files. Each sub list has two entries:
  %
  %   CaseName FileNameTemplate
  %
  CaseList = {
%    { 'RCE50_RECT'      'HDF5/RCE50_RECT/HDF5/<FPREFIX>-a-AS-2012-01-01-000000-g1.h5'      }
%    { 'RCE50_RECT_S300' 'HDF5/RCE50_RECT_S300/HDF5/<FPREFIX>-a-AS-2012-01-01-000000-g1.h5' }
%    { 'RCE50_RECT_S303' 'HDF5/RCE50_RECT_S303/HDF5/<FPREFIX>-a-AS-2012-01-01-000000-g1.h5' }
%    { 'RCE50_SQ'        'HDF5/RCE50_SQ/HDF5/<FPREFIX>-a-AS-2012-01-01-000000-g1.h5'        }
%    { 'RCE50_2D'        'HDF5/RCE50_2D/HDF5/<FPREFIX>-a-AS-2012-01-01-000000-g1.h5'        }
%    { 'RCE_MATT'        'HDF5/RCE_MATT/HDF5/<FPREFIX>-a-AS-2012-01-01-000000-g1.h5'        }
%    { 'RCE_BASE'        'HDF5/RCE_BASE/HDF5/<FPREFIX>-a-AC-2012-01-01-000000-g1.h5'        }
%    { 'RCE_CNTL'        'HDF5/RCE_CNTL/HDF5/<FPREFIX>-RCE_CNTL-AC-2012-01-01-000000-g1.h5'        }
%    { 'RCE_CNTL_NZ75'   'HDF5/RCE_CNTL/HDF5_NZ75/<FPREFIX>-RCE_CNTL-AC-2012-01-01-000000-g1.h5'   }

%    { 'RCE_BASE'        'HDF5/RCE_BASE/HDF5/<FPREFIX>-a-AC-2012-01-01-000000-g1.h5'                  }
%    { 'RCE_EXP_S50LN'   'HDF5/RCE_EXP_S50LN/HDF5/<FPREFIX>-RCE_EXP_S50LN-AC-2012-01-01-000000-g1.h5' }
%    { 'RCE_EXP_S70LY'   'HDF5/RCE_EXP_S70LY/HDF5/<FPREFIX>-RCE_EXP_S70LY-AC-2012-01-01-000000-g1.h5' }
    { 'RCE_EXP_S70LN'   'HDF5/RCE_EXP_S70LN/HDF5/<FPREFIX>-RCE_EXP_S70LN-AC-2012-01-01-000000-g1.h5' }
    { 'RCE_EXP_S70MY'   'HDF5/RCE_EXP_S70MY/HDF5/<FPREFIX>-RCE_EXP_S70MY-AC-2012-01-01-000000-g1.h5' }
   };
  Ncases = length(CaseList);

  SfcLatVname  = 'lat_flux';
  SfcSensVname = 'sens_flux';
  SfcSwdnVname = 'rshort';
  SfcLwdnVname = 'rlong';
  SfcLwupVname = 'rlongup';
  SfcAlbVname  = 'albedt';
  TopSwdnVname = 'swdn';
  TopSwupVname = 'swup';
  TopLwupVname = 'lwup';

  SfcLatFprefix  = 'lat_flux';
  SfcSensFprefix = 'sens_flux';
  SfcSwdnFprefix = 'rshort';
  SfcLwdnFprefix = 'rlong';
  SfcLwupFprefix = 'rlongup';
  SfcAlbFprefix  = 'albedt';
  TopSwdnFprefix = 'swdn_toa';
  TopSwupFprefix = 'swup_toa';
  TopLwupFprefix = 'lwup_toa';

  % output file specs
  OutFprefix = 'eq_meas';
  ThfVname = 'therm_heat_flux';
  QradVname = 'rad_flux_div';

  fprintf('***************************************************************\n');
  fprintf('Generating domain average time series:\n');

  for icase = 1:Ncases
    Case       = CaseList{icase}{1};
    InFileTmpl = CaseList{icase}{2};

    OutFile = sprintf('%s/%s_%s.h5', Ddir, OutFprefix, Case);

    fprintf('  Case: %s\n', Case);
    fprintf('\n');

    % Assemble input file names
    SfcLatFile  = regexprep(InFileTmpl, '<FPREFIX>', SfcLatFprefix);
    SfcSensFile = regexprep(InFileTmpl, '<FPREFIX>', SfcSensFprefix);

    SfcSwdnFile = regexprep(InFileTmpl, '<FPREFIX>', SfcSwdnFprefix);
    SfcLwdnFile = regexprep(InFileTmpl, '<FPREFIX>', SfcLwdnFprefix);
    SfcLwupFile = regexprep(InFileTmpl, '<FPREFIX>', SfcLwupFprefix);
    SfcAlbFile  = regexprep(InFileTmpl, '<FPREFIX>', SfcAlbFprefix);

    TopSwdnFile  = regexprep(InFileTmpl, '<FPREFIX>', TopSwdnFprefix);
    TopSwupFile  = regexprep(InFileTmpl, '<FPREFIX>', TopSwupFprefix);
    TopLwupFile  = regexprep(InFileTmpl, '<FPREFIX>', TopLwupFprefix);

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

    AVG_SFC_LAT  = zeros([ Nt 1 ]);
    AVG_SFC_SENS = zeros([ Nt 1 ]);

    AVG_SFC_SWDN = zeros([ Nt 1 ]);
    AVG_SFC_LWDN = zeros([ Nt 1 ]);
    AVG_SFC_LWUP = zeros([ Nt 1 ]);
    AVG_SFC_SWUP = zeros([ Nt 1 ]);
    AVG_TOP_SWDN = zeros([ Nt 1 ]);
    AVG_TOP_SWUP = zeros([ Nt 1 ]);
    AVG_TOP_LWUP = zeros([ Nt 1 ]);

    AVG_SFC_ALB  = zeros([ Nt 1 ]);

    for it = 1:Nt
      % Read in the variables -> after squeeze will be (y,x)
      % Throw away the lateral boundaries since these are sometimes
      % set to zero
      if (strcmp(Case, 'RCE50_2D'))
        SFC_LAT  = squeeze(SFC_LAT_VAR.data(it,1,2:end-1));
        SFC_SENS = squeeze(SFC_SENS_VAR.data(it,1,2:end-1));
  
        SFC_SWDN = squeeze(SFC_SWDN_VAR.data(it,1,2:end-1));
        SFC_LWDN = squeeze(SFC_LWDN_VAR.data(it,1,2:end-1));
        SFC_LWUP = squeeze(SFC_LWUP_VAR.data(it,1,2:end-1));
        SFC_ALB  = squeeze(SFC_ALB_VAR.data(it,1,2:end-1));
        TOP_SWDN = squeeze(TOP_SWDN_VAR.data(it,1,1,2:end-1));
        TOP_SWUP = squeeze(TOP_SWUP_VAR.data(it,1,1,2:end-1));
        TOP_LWUP = squeeze(TOP_LWUP_VAR.data(it,1,1,2:end-1));
      else
        SFC_LAT  = squeeze(SFC_LAT_VAR.data(it,2:end-1,2:end-1));
        SFC_SENS = squeeze(SFC_SENS_VAR.data(it,2:end-1,2:end-1));
  
        SFC_SWDN = squeeze(SFC_SWDN_VAR.data(it,2:end-1,2:end-1));
        SFC_LWDN = squeeze(SFC_LWDN_VAR.data(it,2:end-1,2:end-1));
        SFC_LWUP = squeeze(SFC_LWUP_VAR.data(it,2:end-1,2:end-1));
        SFC_ALB  = squeeze(SFC_ALB_VAR.data(it,2:end-1,2:end-1));
        TOP_SWDN = squeeze(TOP_SWDN_VAR.data(it,:,2:end-1,2:end-1));
        TOP_SWUP = squeeze(TOP_SWUP_VAR.data(it,:,2:end-1,2:end-1));
        TOP_LWUP = squeeze(TOP_LWUP_VAR.data(it,:,2:end-1,2:end-1));
      end

      % SFC_SWUP = SFC_SWDN * SFC_ALB
      SFC_SWUP = SFC_SWDN .* SFC_ALB;

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
      % Calculate THF, QRAD per column and then take average
      %
      % Save out the averages of all the component quantities, too.
      %
      % Use double precision for the mean calculations since the size of
      % the vectors can be large.
      %
      TempVar = SFC_LAT + SFC_SENS;
      THF(it) = mean(double(TempVar(:)));

      TempVar = ((TOP_SWUP + TOP_LWUP) - (TOP_SWDN)) - ((SFC_SWUP + SFC_LWUP) - (SFC_SWDN + SFC_LWDN));
      QRAD(it) = mean(double(TempVar(:)));

      AVG_SFC_LAT(it) = mean(double(SFC_LAT(:)));
      AVG_SFC_SENS(it) = mean(double(SFC_SENS(:)));

      AVG_SFC_SWDN(it) = mean(double(SFC_SWDN(:)));
      AVG_SFC_LWDN(it) = mean(double(SFC_LWDN(:)));
      AVG_SFC_LWUP(it) = mean(double(SFC_LWUP(:)));
      AVG_SFC_SWUP(it) = mean(double(SFC_SWUP(:)));
      AVG_TOP_SWDN(it) = mean(double(TOP_SWDN(:)));
      AVG_TOP_SWUP(it) = mean(double(TOP_SWUP(:)));
      AVG_TOP_LWUP(it) = mean(double(TOP_LWUP(:)));

      AVG_SFC_ALB(it)  = mean(double(SFC_ALB(:)));

      if (mod(it,10) == 0)
        fprintf ('  Completed timestep %d out of %d\n', it, Nt);
      end
    end
    fprintf('\n');

    fprintf('    Writing: %s (%s, %s)\n', OutFile, ThfVname, QradVname);
    fprintf('\n');

    hdf5write(OutFile, ThfVname, THF); 
    hdf5write(OutFile, QradVname, QRAD, 'WriteMode', 'append'); 

    hdf5write(OutFile, '/avg_sfc_lat',  AVG_SFC_LAT,  'WriteMode', 'append'); 
    hdf5write(OutFile, '/avg_sfc_sens', AVG_SFC_SENS, 'WriteMode', 'append'); 

    hdf5write(OutFile, '/avg_sfc_swdn', AVG_SFC_SWDN, 'WriteMode', 'append'); 
    hdf5write(OutFile, '/avg_sfc_lwdn', AVG_SFC_LWDN, 'WriteMode', 'append'); 
    hdf5write(OutFile, '/avg_sfc_lwup', AVG_SFC_LWUP, 'WriteMode', 'append'); 
    hdf5write(OutFile, '/avg_sfc_swup', AVG_SFC_SWUP, 'WriteMode', 'append'); 
    hdf5write(OutFile, '/avg_top_swdn', AVG_TOP_SWDN, 'WriteMode', 'append'); 
    hdf5write(OutFile, '/avg_top_swup', AVG_TOP_SWUP, 'WriteMode', 'append'); 
    hdf5write(OutFile, '/avg_top_lwup', AVG_TOP_LWUP, 'WriteMode', 'append'); 

    hdf5write(OutFile, '/avg_sfc_alb',  AVG_SFC_ALB,  'WriteMode', 'append'); 

    hdf5write(OutFile, '/t_coords', T,      'WriteMode', 'append');
  end
end
