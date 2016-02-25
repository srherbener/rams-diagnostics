function [ ] = GenTotalMassMeas()
% GenTotalMassMeas function to calculate total mass in volume or area

  % Cases
  CaseList = {
    'TSD_SAL_DUST'
%    'TSD_SAL_NODUST'
%    'TSD_NONSAL_DUST'
%    'TSD_NONSAL_NODUST'
    };
  Ncases = length(CaseList);

  % Description of measurements
  % third argument is scaling factor to convert to g/m2 or g/m3
  FileList = {
    { 'DIAGS/hda_meas_ts_dust_<CASE>.h5'   '/sal_sum_dust_sfc'      1e-2  '/sal_dust_sfc_total_mass'     5500 }  % units are ug/cm2 so scale by 1e-2 to get g/m2
    { 'DIAGS/hda_meas_ts_dust_<CASE>.h5'   '/sal_sum_dust_mass'     1e-6  '/sal_dust_total_mass'         5500 }  % units are ug/m3 so scale by 1e-6 to get g/m3
    { 'DIAGS/hda_meas_ts_dust_<CASE>.h5'   '/sal_sum_tracer_mass'   1e-6  '/sal_tracer_total_mass'       5500 }  % units are ug/m3 so scale by 1e-6 to get g/m3
    { 'DIAGS/hda_meas_ts_dust_<CASE>.h5'   '/sal_sum_dust_hydro'    1e-6  '/sal_dust_hydro_total_mass'   5500 }  % units are ug/m3 so scale by 1e-6 to get g/m3
    { 'DIAGS/hda_meas_ts_aero_<CASE>.h5'   '/sal_sum_aero_mass'     1e-6  '/sal_aero_total_mass'         5500 }  % units are ug/m3 so scale by 1e-6 to get g/m3
    { 'DIAGS/hda_meas_ts_ccn_<CASE>.h5'    '/sal_sum_ccn_mass'      1e-6  '/sal_ccn_total_mass'          5500 }  % units are ug/m3 so scale by 1e-6 to get g/m3
    { 'DIAGS/hda_meas_ts_ra_<CASE>.h5'     '/sal_sum_ra_mass'       1e-6  '/sal_ra_total_mass'           5500 }  % units are ug/m3 so scale by 1e-6 to get g/m3

    { 'DIAGS/hda_meas_ts_cond_<CASE>.h5'   '/sal_sum_tcond_mass'    'rho' '/sal_tcond_total_mass'        5500 }  % units are g/kg so scale by rho (density) to get g/m3
    { 'DIAGS/hda_meas_ts_precip_<CASE>.h5' '/sal_sum_accpcp_mass'   1e3   '/sal_accpcp_total_mass'       5500 }  % units are mm (=kg/m2, assuming density of water
                                                                                                                 %   to be 1000 kg/m3) so scale by 1e3 to get g/m2
 
    { 'DIAGS/hda_meas_ts_dust_<CASE>.h5'   '/spath_sum_dust_sfc'    1e-2 '/spath_dust_sfc_total_mass'   5500 }  % units are ug/cm2 so scale by 1e-2 to get g/m2
    { 'DIAGS/hda_meas_ts_dust_<CASE>.h5'   '/spath_sum_dust_mass'   1e-6 '/spath_dust_total_mass'       5500 }  % units are ug/m3 so scale by 1e-6 to get g/m3
    { 'DIAGS/hda_meas_ts_dust_<CASE>.h5'   '/spath_sum_tracer_mass' 1e-6 '/spath_tracer_total_mass'     5500 }  % units are ug/m3 so scale by 1e-6 to get g/m3
    { 'DIAGS/hda_meas_ts_dust_<CASE>.h5'   '/spath_sum_dust_hydro'  1e-6 '/spath_dust_hydro_total_mass' 5500 }  % units are ug/m3 so scale by 1e-6 to get g/m3
    { 'DIAGS/hda_meas_ts_aero_<CASE>.h5'   '/spath_sum_aero_mass'   1e-6 '/spath_aero_total_mass'       5500 }  % units are ug/m3 so scale by 1e-6 to get g/m3
    { 'DIAGS/hda_meas_ts_ccn_<CASE>.h5'    '/spath_sum_ccn_mass'    1e-6 '/spath_ccn_total_mass'        5500 }  % units are ug/m3 so scale by 1e-6 to get g/m3
    { 'DIAGS/hda_meas_ts_ra_<CASE>.h5'     '/spath_sum_ra_mass'     1e-6 '/spath_ra_total_mass'         5500 }  % units are ug/m3 so scale by 1e-6 to get g/m3

    { 'DIAGS/hist_sums_az_dust_<CASE>.h5'  '/all_sum_dust_sfc'      1e-2 '/storm_dust_sfc_total_mass'   5500 }  % units are ug/cm2 so scale by 1e-2 to get g/m2
    { 'DIAGS/hist_sums_az_dust_<CASE>.h5'  '/all_sum_dust_mass'     1e-6 '/storm_dust_total_mass'       5500 }  % units are ug/m3 so scale by 1e-6 to get g/m3
    { 'DIAGS/hist_sums_az_dust_<CASE>.h5'  '/all_sum_tracer_mass'   1e-6 '/storm_tracer_total_mass'     5500 }  % units are ug/m3 so scale by 1e-6 to get g/m3
    { 'DIAGS/hist_sums_az_dust_<CASE>.h5'  '/all_sum_dust_hydro'    1e-6 '/storm_dust_hydro_total_mass' 5500 }  % units are ug/m3 so scale by 1e-6 to get g/m3
    { 'DIAGS/hist_sums_az_aero_<CASE>.h5'  '/all_sum_aero_mass'     1e-6 '/storm_aero_total_mass'       5500 }  % units are ug/m3 so scale by 1e-6 to get g/m3
    { 'DIAGS/hist_sums_az_ccn_<CASE>.h5'   '/all_sum_ccn_mass'      1e-6 '/storm_ccn_total_mass'        5500 }  % units are ug/m3 so scale by 1e-6 to get g/m3
    { 'DIAGS/hist_sums_az_ra_<CASE>.h5'    '/all_sum_ra_mass'       1e-6 '/storm_ra_total_mass'         5500 }  % units are ug/m3 so scale by 1e-6 to get g/m3
    };
  Nfiles = length(FileList);

  HorizArea = 9000 * 9000; % horizontal grid cell area, 3km X 3km

  % Air density - RAMS reference state (DN101D) from the TS Debby runs, kg/m3
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

  for icase = 1:Ncases
    Case = CaseList{icase};

    fprintf('*****************************************************************\n');
    fprintf('Generating total mass measurements:\n');
    fprintf('  Case: %s\n', Case);
    fprintf('\n');

    % Place all measurements into one case specific output file.
    % If the file exists, remove it so that the HDF5 commands
    % can effectively re-create datasets.
    OutFile = sprintf('DIAGS/total_mass_%s.h5', Case);
    if (exist(OutFile, 'file') == 2)
      delete(OutFile);
    end

    for ifile = 1:Nfiles
      Ifile    = FileList{ifile}{1};
      Ivname   = FileList{ifile}{2};
      Scale    = FileList{ifile}{3};
      OutVname = FileList{ifile}{4};
      Zlev     = FileList{ifile}{5};

      InFile  = regexprep(Ifile, '<CASE>', Case);

      fprintf('  Reading: %s (%s)\n', InFile, Ivname);

      % HDATA will be (x,y,t)
      HDATA = squeeze(h5read(InFile, Ivname));
      X = squeeze(h5read(InFile, '/x_coords'));
      Y = squeeze(h5read(InFile, '/y_coords'));
      Z = squeeze(h5read(InFile, '/z_coords'));
      T = squeeze(h5read(InFile, '/t_coords'));

      Nx = length(X);
      Ny = length(Y);
      Nz = length(Z);
      Nt = length(T);

      Hsize = size(HDATA);

      % Apply the scaling factor
      if (strcmp(Scale, 'rho'))
        if((Hsize(1) > 1) && (Hsize(2) > 1))
          % Hdata is (z,t)
          HDATA = HDATA .* repmat(RhoAir, [ 1 Nt ]);
        else
          % Hdata is (t);
          HDATA = HDATA .* RhoAir(2);
        end
      else
        HDATA = HDATA .* Scale;
      end

      % HDATA is either (z,t) or (t)
      Z1 = find(Z >= Zlev, 1, 'first');
      if ((Hsize(1) > 1) && (Hsize(2) > 1))
        % HDATA is (z,t)
        % Create a single column with delta z values
        %   Just repeat delta z for top grid cell (should be
        %   the same given RAMS max delta z configuration)
        DELTA_Z = Z(2:end) - Z(1:end-1);
        DELTA_Z(Nz) = DELTA_Z(Nz-1);

        % Repeat the column Nt times to create delta z matching
        % up with the size of HDATA.
        DELTA_Z = repmat(DELTA_Z, [ 1 Nt ]);

        % Do the vertical integration (summation)
        MEAS = sum(HDATA .* DELTA_Z);
        MEAS_LLEV = sum(HDATA(1:Z1-1,:) .* DELTA_Z(1:Z1-1,:));
        MEAS_HLEV = sum(HDATA(Z1:end,:) .* DELTA_Z(Z1:end,:));
      else
        MEAS = HDATA;
        MEAS_LLEV = HDATA; % dummy placeholder since no z structure
        MEAS_HLEV = HDATA; % dummy placeholder since no z structure
      end

      % Do the horizontal integration
      MEAS = MEAS .* HorizArea;
      MEAS_LLEV = MEAS_LLEV .* HorizArea;
      MEAS_HLEV = MEAS_HLEV .* HorizArea;

      % If this is the first set, write out the coordinates into the output file
      % so that subsequent variables can have these coordinates attached.
      if (ifile == 1)
        Xname = '/x_coords';
        Yname = '/y_coords';
        Zname = '/z_coords';
        Tname = '/t_coords';
  
        CreateDimensionsXyzt(OutFile, X, Y, Z, T, Xname, Yname, Zname, Tname);
        % Add COARDS annotations
        NotateDimensionsXyzt(OutFile, Xname, Yname, Zname, Tname);
      end

      Vsize = Nt;
      DimOrder = { 't' };

      fprintf('  Writing: %s (%s)\n', OutFile, OutVname);
      h5create(OutFile, OutVname, Vsize);
      h5write(OutFile, OutVname, MEAS);

      OutLlevVname = sprintf('%s_llev', OutVname);
      fprintf('  Writing: %s (%s)\n', OutFile, OutLlevVname);
      h5create(OutFile, OutLlevVname, Vsize);
      h5write(OutFile, OutLlevVname, MEAS_LLEV);

      OutHlevVname = sprintf('%s_hlev', OutVname);
      fprintf('  Writing: %s (%s)\n', OutFile, OutHlevVname);
      h5create(OutFile, OutHlevVname, Vsize);
      h5write(OutFile, OutHlevVname, MEAS_HLEV);

      AttachDimensionsXyzt(OutFile, OutVname, DimOrder, Xname, Yname, Zname, Tname);
      AttachDimensionsXyzt(OutFile, OutLlevVname, DimOrder, Xname, Yname, Zname, Tname);
      AttachDimensionsXyzt(OutFile, OutHlevVname, DimOrder, Xname, Yname, Zname, Tname);

      fprintf('\n');
    end
  end
end 
