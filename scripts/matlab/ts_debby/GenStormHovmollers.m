function [ ] = GenStormHovmollers()

  Ddir = 'DIAGS';
  if (exist(Ddir, 'dir') ~= 7)
    mkdir(Ddir);
  end

  Cases = {
   'TSD_SD'
   'TSD_SD_1G'
%   'TSD_SAL_DUST'
%   'TSD_SAL_NODUST'
%   'TSD_NONSAL_DUST'
%   'TSD_NONSAL_NODUST'
   };
  Ncases = length(Cases);

  % From TSD_SAL_DUST RAMS output (Reference density)
  RhoAir = [
    1.196
    1.191
    1.185
    1.179
    1.172
    1.165
    1.157
    1.147
    1.136
    1.124
    1.111
    1.097
    1.082
    1.066
    1.051
    1.034
    1.017
    1.000
    0.982
    0.963
    0.945
    0.925
    0.905
    0.883
    0.860
    0.833
    0.804
    0.773
    0.741
    0.711
    0.679
    0.647
    0.612
    0.577
    0.541
    0.505
    0.469
    0.432
    0.394
    0.355
    0.316
    0.279
    0.244
    0.210
    0.179
    0.150
    0.126
    0.105
    0.087
    0.073
    0.062
    0.052
    0.044
    0.038
    0.032
    0.027
    ];

  % Description of hovmollers
  HovmollerList = {
    % in_file in_var out_var pre_sal_rb_var sal_core_var sal_rb_var

%    % Region: All
%    { 'DIAGS/hist_meas_az_aero_<CASE>.h5' '/all_whole_aero_mass_ts' '/all_whole_aero_mass' 0.01 1000 }

%    % Region: SPATH
%    { 'DIAGS/hist_meas_ts_dust_hydro_<CASE>.h5' '/spath_dust_cloud_ts' '/spath_dust_cloud' 0.001 1000 }
%    { 'DIAGS/hist_meas_ts_dust_hydro_<CASE>.h5' '/spath_dust_rain_ts'  '/spath_dust_rain'  0.001 1000 }
%    { 'DIAGS/hist_meas_ts_dust_hydro_<CASE>.h5' '/spath_dust_pris_ts'  '/spath_dust_pris'  0.001 1000 }
%    { 'DIAGS/hist_meas_ts_dust_hydro_<CASE>.h5' '/spath_dust_snow_ts'  '/spath_dust_snow'  0.001 1000 }
%    { 'DIAGS/hist_meas_ts_dust_hydro_<CASE>.h5' '/spath_dust_aggr_ts'  '/spath_dust_aggr'  0.001 1000 }
%    { 'DIAGS/hist_meas_ts_dust_hydro_<CASE>.h5' '/spath_dust_graup_ts' '/spath_dust_graup' 0.001 1000 }
%    { 'DIAGS/hist_meas_ts_dust_hydro_<CASE>.h5' '/spath_dust_hail_ts'  '/spath_dust_hail'  0.001 1000 }
%
%    { 'DIAGS/hist_meas_ts_dust_<CASE>.h5'       '/spath_d1_mass_ts'     '/spath_d1_mass'     0.001 1000 }
%    { 'DIAGS/hist_meas_ts_dust_<CASE>.h5'       '/spath_d2_mass_ts'     '/spath_d2_mass'     0.01  1000 }
%    { 'DIAGS/hist_meas_ts_dust_<CASE>.h5'       '/spath_dust_mass_ts'   '/spath_dust_mass'   0.01  1000 }
%    { 'DIAGS/hist_meas_ts_dust_<CASE>.h5'       '/spath_tracer_mass_ts' '/spath_tracer_mass' 0.01  1000 }
%    { 'DIAGS/hist_meas_ts_dust_hydro_<CASE>.h5' '/spath_dust_hydro_ts'  '/spath_dust_hydro'  0.001 1000 }
%    { 'DIAGS/hist_meas_ts_ccn_<CASE>.h5'        '/spath_ccn_mass_ts'    '/spath_ccn_mass'    0.001 1000 }
%    { 'DIAGS/hist_meas_ts_ra_<CASE>.h5'         '/spath_ra_mass_ts'     '/spath_ra_mass'     0.01  1000 }
%    { 'DIAGS/hist_meas_ts_aero_<CASE>.h5'       '/spath_aero_mass_ts'   '/spath_aero_mass'   0.01  1000 }
%
%    { 'DIAGS/hist_meas_ts_dust_<CASE>.h5'       '/sal_d1_mass_ts'     '/sal_d1_mass'     0.001 1000 }
%    { 'DIAGS/hist_meas_ts_dust_<CASE>.h5'       '/sal_d2_mass_ts'     '/sal_d2_mass'     0.01  1000 }
%    { 'DIAGS/hist_meas_ts_dust_<CASE>.h5'       '/sal_dust_mass_ts'   '/sal_dust_mass'   0.01  1000 }
%    { 'DIAGS/hist_meas_ts_dust_<CASE>.h5'       '/sal_tracer_mass_ts' '/sal_tracer_mass' 0.01  1000 }
%    { 'DIAGS/hist_meas_ts_dust_hydro_<CASE>.h5' '/sal_dust_hydro_ts'  '/sal_dust_hydro'  0.001 1000 }
%    { 'DIAGS/hist_meas_ts_ra_<CASE>.h5'         '/sal_ra_mass_ts'     '/sal_ra_mass'     0.01  1000 }
%    { 'DIAGS/hist_meas_ts_ccn_<CASE>.h5'        '/sal_ccn_mass_ts'    '/sal_ccn_mass'    0.001 1000 }
%    { 'DIAGS/hist_meas_ts_aero_<CASE>.h5'       '/sal_aero_mass_ts'   '/sal_aero_mass'   0.01  1000 }

%    { 'DIAGS/hist_meas_ts_dust_<CASE>.h5'       '/sal_ar_d1_mass_ts'     '/sal_ar_d1_mass'     0.001 1000 }
%    { 'DIAGS/hist_meas_ts_dust_<CASE>.h5'       '/sal_ar_d2_mass_ts'     '/sal_ar_d2_mass'     0.01  1000 }
    { 'DIAGS/hist_meas_ts_dust_<CASE>.h5'       '/sal_ar_dust_mass_ts'   '/sal_ar_dust_mass'   0.01  1000 }
%    { 'DIAGS/hist_meas_ts_dust_<CASE>.h5'       '/sal_ar_tracer_mass_ts' '/sal_ar_tracer_mass' 0.01  1000 }
    { 'DIAGS/hist_meas_ts_dust_hydro_<CASE>.h5' '/sal_ar_dust_hydro_ts'  '/sal_ar_dust_hydro'  0.001 1000 }
    { 'DIAGS/hist_meas_ts_ra_<CASE>.h5'         '/sal_ar_ra_mass_ts'     '/sal_ar_ra_mass'     0.01  1000 }
%    { 'DIAGS/hist_meas_ts_ccn_<CASE>.h5'        '/sal_ar_ccn_mass_ts'    '/sal_ar_ccn_mass'    0.001 1000 }
%    { 'DIAGS/hist_meas_ts_aero_<CASE>.h5'       '/sal_ar_aero_mass_ts'   '/sal_ar_aero_mass'   0.01  1000 }

%    { 'DIAGS/hist_meas_az_dust_<CASE>.h5'       '/all_whole_d1_mass_ts'     '/storm_d1_mass'     0.001 1000 }
%    { 'DIAGS/hist_meas_az_dust_<CASE>.h5'       '/all_whole_d2_mass_ts'     '/storm_d2_mass'     0.01  1000 }
%    { 'DIAGS/hist_meas_az_dust_<CASE>.h5'       '/all_whole_dust_mass_ts'   '/storm_dust_mass'   0.01  1000 }
%    { 'DIAGS/hist_meas_az_dust_<CASE>.h5'       '/all_whole_tracer_mass_ts' '/storm_tracer_mass' 0.01  1000 }
%    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_whole_dust_hydro_ts'  '/storm_dust_hydro'  0.001 1000 }
%    { 'DIAGS/hist_meas_az_ra_<CASE>.h5'         '/all_whole_ra_mass_ts'     '/storm_ra_mass'     0.01  1000 }
%    { 'DIAGS/hist_meas_az_ccn_<CASE>.h5'        '/all_whole_ccn_mass_ts'    '/storm_ccn_mass'    0.001 1000 }
%    { 'DIAGS/hist_meas_az_aero_<CASE>.h5'       '/all_whole_aero_mass_ts'   '/storm_aero_mass'   0.01  1000 }
  };
  Nsets = length(HovmollerList);


  for icase = 1:Ncases
    Case = Cases{icase};

    fprintf('**********************************************************\n');
    fprintf('Generating storm hovmollers for case: %s\n', Case);
    fprintf('\n');

    % Place all vars into one case specific output file.
    % If the file exists, remove it so that the HDF5 commands
    % can effectively re-create datasets.
    OutFile = sprintf('DIAGS/storm_hovmollers_%s.h5', Case);
    if (exist(OutFile, 'file') == 2)
      delete(OutFile);
    end

    icount = 0;
    for iset = 1:Nsets
      InFile   = regexprep(HovmollerList{iset}{1}, '<CASE>', Case);
      InVname  = HovmollerList{iset}{2};
      OutVname = HovmollerList{iset}{3};
      Zmin     = HovmollerList{iset}{4};
      Zmax     = HovmollerList{iset}{5};

      OnFinalProf = ~isempty(regexp(InVname, '_f_'));

      % skip this hovmoller set if doing dust and on a NODUST case
      if ((~isempty(regexp(Case, 'NODUST'))) && ...
          ((~isempty(regexp(InVname, 'd[12]_num'))) || ...
           (~isempty(regexp(InVname, 'd[12]_mass'))) || ...
           (~isempty(regexp(InVname, '_dust_'))) || ...
           (~isempty(regexp(InVname, '_dustifn_'))) || ...
           (~isempty(regexp(InVname, '_trdust[12]_'))) || ...
           (~isempty(regexp(InVname, '_tracer[12]')))))
        continue
      else
        icount = icount + 1;
      end

      % Read in the variables, split into pre-SAL and SAL sections
      fprintf('  Reading: %s (%s)\n', InFile, InVname);

      VAR = squeeze(h5read(InFile, InVname));
      X   = squeeze(h5read(InFile, '/x_coords'));
      Y   = squeeze(h5read(InFile, '/y_coords'));
      Z   = squeeze(h5read(InFile, '/z_coords'));
      T   = squeeze(h5read(InFile, '/t_coords'));

      Nx = length(X);
      Ny = length(Y);
      Nz = length(Z);
      Nt = length(T);

      % If a hydrometeor number concentration, then convert #/kg to #/cm^2.
      if (~isempty(regexp(InVname, 'cloud_num')) || ...
          ~isempty(regexp(InVname, 'rain_num'))  || ...
          ~isempty(regexp(InVname, 'pris_num'))  || ...
          ~isempty(regexp(InVname, 'snow_num'))  || ...
          ~isempty(regexp(InVname, 'aggr_num'))  || ...
          ~isempty(regexp(InVname, 'graup_num')) || ...
          ~isempty(regexp(InVname, 'hail_num')))
        VAR      = VAR .* RhoAir .* 1e-6;
      end

      % Change values outside the range (Zmin,Zmax) to nan's to help make the contour
      % plots more readable
      VAR(VAR < Zmin) = nan;
      VAR(VAR > Zmax) = nan;

      % If this is the first set, write out the coordinates into the output file
      % so that subsequent variables can have these coordinates attached.
      if (icount == 1)
        Xname = '/x_coords';
        Yname = '/y_coords';
        Zname = '/z_coords';
        Tname = '/t_coords';
  
        CreateDimensionsXyzt(OutFile, X, Y, Z, T, Xname, Yname, Zname, Tname);
        % Add COARDS annotations
        NotateDimensionsXyzt(OutFile, Xname, Yname, Zname, Tname);
      end

      Vsize = size(VAR);
      DimOrder = { 'z' 't' };

      % Collect output vars into one case specific file.
      fprintf('  Writing: %s (%s)\n', OutFile, OutVname);
      h5create(OutFile, OutVname, Vsize);
      h5write(OutFile, OutVname, VAR);
      AttachDimensionsXyzt(OutFile, OutVname, DimOrder, Xname, Yname, Zname, Tname);

      fprintf('\n');
    end
  end
end
