function [ ] = GenStormMeas()
% GenStormMeas function to create storm measurements (max wind, min pressure, rmw, ike, etc)

  Ddir = 'DIAGS';

  % make sure output directory exists
  if (exist(Ddir, 'dir') ~= 7)
    mkdir(Ddir);
  end

  % Temporal ranges
  TstartPreSal = 10;
  TendPreSal   = 30;

  TstartSal = 40;
  TendSal   = 60;

  % Cases
  CaseList = {
    'TSD_SAL_DUST'
    'TSD_SAL_NODUST'
    'TSD_NONSAL_DUST'
    'TSD_NONSAL_NODUST'
    };

  Ncases = length(CaseList);

  % Description of measurements
  MeasList = {
    % vertical coord is height
    { 'min_slp'    'DIAGS/hist_meas_az_press_<CASE>.h5' '/all_press_sfc_ts'    0 250 0 1000 'min' }
    { 'ps_min_slp' 'DIAGS/hist_meas_az_press_<CASE>.h5' '/all_ps_press_sfc_ts' 0 250 0 1000 'min' }
    { 's_min_slp'  'DIAGS/hist_meas_az_press_<CASE>.h5' '/all_s_press_sfc_ts'  0 250 0 1000 'min' }

    { 'max_wind_t'    'DIAGS/hist_meas_az_speed_<CASE>.h5' '/all_speed_t_maxlev_ts'    0 250 0 2000 'max' }
    { 'ps_max_wind_t' 'DIAGS/hist_meas_az_speed_<CASE>.h5' '/all_ps_speed_t_maxlev_ts' 0 250 0 2000 'max' }
    { 's_max_wind_t'  'DIAGS/hist_meas_az_speed_<CASE>.h5' '/all_s_speed_t_maxlev_ts'  0 250 0 2000 'max' }

    { 'avg_wind_t'    'DIAGS/hist_meas_az_speed_<CASE>.h5' '/all_speed_t_ts'    0 250 0 2000 'avg' }
    { 'ps_avg_wind_t' 'DIAGS/hist_meas_az_speed_<CASE>.h5' '/all_ps_speed_t_ts' 0 250 0 2000 'avg' }
    { 's_avg_wind_t'  'DIAGS/hist_meas_az_speed_<CASE>.h5' '/all_s_speed_t_ts'  0 250 0 2000 'avg' }

    { 'rmw_t'       'DIAGS/hist_meas_az_speed_<CASE>.h5'   '/all_speed_t_maxlev_ts'    0 250 0 2000 'rmw' }
    { 'ps_rmw_t'    'DIAGS/hist_meas_az_speed_<CASE>.h5'   '/all_ps_speed_t_maxlev_ts' 0 250 0 2000 'rmw' }
    { 's_rmw_t'     'DIAGS/hist_meas_az_speed_<CASE>.h5'   '/all_s_speed_t_maxlev_ts'  0 250 0 2000 'rmw' }

    { 'ike'      'TsAveragedData/horiz_ke_<CASE>.h5' '/horiz_ke'    0 250 0 2000 'na'  }
    { 'ps_ike'   'TsAveragedData/horiz_ke_<CASE>.h5' '/horiz_ke'    0 250 0 2000 'na'  }
    { 's_ike'    'TsAveragedData/horiz_ke_<CASE>.h5' '/horiz_ke'    0 250 0 2000 'na'  }

    { 'pcprate'     'DIAGS/hist_meas_az_pcprate_<CASE>.h5' '/all_pcprate_ts'     0 250 0 1000 'max' }
    { 'ps_pcprate'  'DIAGS/hist_meas_az_pcprate_<CASE>.h5' '/all_ps_pcprate_ts'  0 250 0 1000 'max' }
    { 's_pcprate'   'DIAGS/hist_meas_az_pcprate_<CASE>.h5' '/all_s_pcprate_ts'   0 250 0 1000 'max' }

    % vertical coord is pressure
    { 'max_wind_t_p'      'DIAGS/hist_meas_az_p_speed_<CASE>.h5' '/all_speed_t_maxlev_ts'    0 250 1010 800 'max' }
    { 'ps_max_wind_t_p'   'DIAGS/hist_meas_az_p_speed_<CASE>.h5' '/all_ps_speed_t_maxlev_ts' 0 250 1010 800 'max' }
    { 's_max_wind_t_p'    'DIAGS/hist_meas_az_p_speed_<CASE>.h5' '/all_s_speed_t_maxlev_ts'  0 250 1010 800 'max' }

    { 'avg_wind_t_p'      'DIAGS/hist_meas_az_p_speed_<CASE>.h5' '/all_speed_t_ts'    0 250 1010 800 'avg' }
    { 'ps_avg_wind_t_p'   'DIAGS/hist_meas_az_p_speed_<CASE>.h5' '/all_ps_speed_t_ts' 0 250 1010 800 'avg' }
    { 's_avg_wind_t_p'    'DIAGS/hist_meas_az_p_speed_<CASE>.h5' '/all_s_speed_t_ts'  0 250 1010 800 'avg' }

    { 'rmw_t_p'         'DIAGS/hist_meas_az_p_speed_<CASE>.h5'   '/all_speed_t_maxlev_ts'    0 250 1010 800 'rmw' }
    { 'ps_rmw_t_p'      'DIAGS/hist_meas_az_p_speed_<CASE>.h5'   '/all_ps_speed_t_maxlev_ts' 0 250 1010 800 'rmw' }
    { 's_rmw_t_p'       'DIAGS/hist_meas_az_p_speed_<CASE>.h5'   '/all_s_speed_t_maxlev_ts'  0 250 1010 800 'rmw' }

    };

  for icase = 1:Ncases
    Case = CaseList{icase};

    fprintf('*****************************************************************\n');
    fprintf('Generating storm measurements:\n');
    fprintf('  Case: %s\n', Case);
    fprintf('\n');

    % Write out measurement
    OutFile = sprintf('%s/storm_meas_%s.h5', Ddir, Case);
    fprintf('    Writing: %s\n', OutFile)
    fprintf('\n');

    % get old file out of the way
    if (exist(OutFile, 'file') == 2)
      delete(OutFile);
    end

    for imeas = 1:length(MeasList)
      Mname    = MeasList{imeas}{1};
      Mfile    = MeasList{imeas}{2};
      Mvname   = MeasList{imeas}{3};
      Rmin     = MeasList{imeas}{4} .* 1000;  % convert km to m
      Rmax     = MeasList{imeas}{5} .* 1000;  % convert km to m
      Zbot     = MeasList{imeas}{6};
      Ztop     = MeasList{imeas}{7};
      Mmeas    = MeasList{imeas}{8};

      Mfile = regexprep(Mfile, '<CASE>', Case);

      fprintf('  Measurement: %s\n', Mname);
      fprintf('    Reading: %s (%s)\n', Mfile, Mvname);
      fprintf('    Rmin: %.2f\n', Rmin);
      fprintf('    Rmax: %.2f\n', Rmax);
      fprintf('    Zbot: %.2f\n', Zbot);
      fprintf('    Ztop: %.2f\n', Ztop);
      fprintf('    Measurement: %s\n', Mmeas);
      fprintf('\n');

      % MDATA will be either 1D, 2D or 3D
      %    1D: (t)
      %    2D: (r,t)
      %    3D: (r,z,t)
      MDATA = squeeze(h5read(Mfile, Mvname));
      R = squeeze(h5read(Mfile, '/x_coords'));
      Y = squeeze(h5read(Mfile, '/y_coords'));
      Z = squeeze(h5read(Mfile, '/z_coords'));
      T = squeeze(h5read(Mfile, '/t_coords'));

      Nx = length(R);
      Ny = length(Y);
      Nz = length(Z);
      Nt = length(T);

      R1 = find(R >= Rmin, 1, 'first');
      R2 = find(R <= Rmax, 1, 'last');

      % Crude test to distinguish pressure and height for vertical
      % coordinates. Pressure values decrease from bottom to top;
      % height values increase.
      Z1 = find(Z >= Zbot, 1, 'first');
      Z2 = find(Z <= Ztop, 1, 'last');
      if (Nz > 1)
        if (Z(1) > Z(2))
          Z1 = find(Z <= Zbot, 1, 'first');
          Z2 = find(Z >= Ztop, 1, 'last');
        end
      end

      % If first measurement, write out coordinates for subsequent
      % variable attaching.
      if (imeas == 1)
        Xname = '/x_coords';
        Yname = '/y_coords';
        Zname = '/z_coords';
        Tname = '/t_coords';

        CreateDimensionsXyzt(OutFile, R, Y, Z, T, Xname, Yname, Zname, Tname);
        % Add COARDS annotations
        NotateDimensionsXyzt(OutFile, Xname, Yname, Zname, Tname);
      end

      % If doing ps_ike or s_ike, select out the pre-SAL or SAL time periods
      switch(Mname)
        case { 'ps_ike' 's_ike' }
          ST = (T ./ 3600) - 42; % simulation time
          if (strcmp(Mname, 'ps_ike'))
            T1 = find(ST >= TstartPreSal, 1, 'first');
            T2 = find(ST <= TendPreSal,   1, 'last');
          else
            T1 = find(ST >= TstartSal, 1, 'first');
            T2 = find(ST <= TendSal,   1, 'last');
          end

          MDATA = MDATA(T1:T2);
      end

      % Grab the number of dimensions. When you have a vector (size is [ 1 n ] or [ n 1 ]),
      % Ndims will be set to 2 when you really want it to be 1.
      Ndims = ndims(MDATA);
      if (Ndims == 2)
        if ((size(MDATA,1) == 1) || (size(MDATA,2) == 1))
          Ndims = 1;
        end
      end

      % Some inputs may have a shorter range on dimensions
      % compared to the full domain and time period
      if (Ndims == 1)
        VarNt = length(MDATA);
      else
        TempSize = size(MDATA);
        VarNt = TempSize(end);
      end

      % Reduce data to 1D time series using Mmeas.
      %  If data is 1D, then done
      %  If data is 2D, run Mmeas on 1st dimension
      %  If data is 3D, run Mmeas on 1st, 2nd dimensions
      if (Ndims > 1)
        Tseries = zeros([ 1 VarNt ]);
        for it = 1:VarNt
          if (ndims(MDATA) == 2)
             % 2D -> (r,t)
             %   M_TEST <- (r)
             M_TEST = squeeze(MDATA(R1:R2,it));
          elseif (ndims(MDATA) == 3)
             % 3D -> (r,z,t)
             %   M_TEST <- (r,z)
             M_TEST = squeeze(MDATA(R1:R2,Z1:Z2,it));
          end

          switch(Mmeas)
            case { 'min' }
              Tseries(it) = min(M_TEST(:));

            case { 'max' 'rmw' }
              % record the index and value of the maximum
              % if doing rmw, use the index to look up the radius
              %
              [ Val Ind ] = max(M_TEST(:));
  
              if (strcmp(Mmeas, 'max'))
                Tseries(it) = Val(1);
              elseif (strcmp(Mmeas, 'rmw'))
                [ Rind Zind ] = ind2sub(size(M_TEST), Ind(1));
                Tseries(it) = R(Rind);
              end

            case { 'avg' }
              Tseries(it) = mean(M_TEST(:));

          end
        end
      else
        Tseries = MDATA;
      end

      % Smooth out RMW time series
      if (strcmp(Mmeas, 'rmw'))
        Tseries = SmoothFillTseries(Tseries, VarNt, 5);
      end

      % Write out each measurement (time series)
      Ovname = sprintf('/%s', Mname);
      DimOrder = { 't' };  % so far, all measurement are time series
      Vsize = VarNt;

      h5create(OutFile, Ovname, Vsize);
      h5write (OutFile, Ovname, Tseries);

      % Attach dimensions
      AttachDimensionsXyzt(OutFile, Ovname, DimOrder, Xname, Yname, Zname, Tname);
    end
  end
end 
