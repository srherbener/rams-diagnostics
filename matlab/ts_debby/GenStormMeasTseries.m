function [ ] = GenStormMeasTseries()
% GenStormMeasTseries function to create time series of storm measurements

  Ddir = 'DIAGS';

  % make sure output directory exists
  if (exist(Ddir, 'dir') ~= 7)
    mkdir(Ddir);
  end

  % Cases
  CaseList = {
    'TSD_SAL_DUST'
    'TSD_SAL_NODUST'
    'TSD_NONSAL_DUST'
    'TSD_NONSAL_NODUST'
    };

  Ncases = length(CaseList);

  % Select points that fall between 0 and 250 km radii
  R1 = 1;
  R2 = 28;

  % Select points that are betwen 0 and 1000 m height
  Z1 = 2;
  Z2 = 14;

  % Description of measurements
  MeasList = {
    { 'min_slp_fa'  'DIAGS/hist_meas_press_<CASE>.h5' '/avg_sea_press_fa' 2 '(R1:R2,:)'      'min' }
    { 'min_slp_wm'  'DIAGS/hist_meas_press_<CASE>.h5' '/avg_sea_press_wm' 2 '(R1:R2,:)'      'min' }

    { 'max_wind_fa'     'DIAGS/hist_meas_speed_<CASE>.h5' '/avg_speed_t_fa'  3  '(R1:R2,Z1:Z2,:)' 'max' }
    { 'max_wind_wm'     'DIAGS/hist_meas_speed_<CASE>.h5' '/avg_speed_t_wm'  3  '(R1:R2,Z1:Z2,:)' 'max' }
    { 'max_wind_10m_fa' 'DIAGS/hist_meas_speed_<CASE>.h5' '/avg_speed10m_fa' 2  '(R1:R2,:)'      'max' }
    { 'max_wind_10m_wm' 'DIAGS/hist_meas_speed_<CASE>.h5' '/avg_speed10m_wm' 2  '(R1:R2,:)'      'max' }

    { 'ike'             'TsAveragedData/horiz_ke_<CASE>.h5' '/horiz_ke' 1  '(:)'      'na' }

    { 'rmw_fa'          'DIAGS/hist_meas_speed_<CASE>.h5'  '/avg_speed_t_fa'  3 '(R1:R2,Z1:Z2,:)'  'rmw' }
    { 'rmw_wm'          'DIAGS/hist_meas_speed_<CASE>.h5'  '/avg_speed_t_wm'  3 '(R1:R2,Z1:Z2,:)'  'rmw' }
    { 'rmw_10m_fa'      'DIAGS/hist_meas_speed_<CASE>.h5'  '/avg_speed10m_fa' 2 '(R1:R2,:)'        'rmw' }
    { 'rmw_10m_wm'      'DIAGS/hist_meas_speed_<CASE>.h5'  '/avg_speed10m_wm' 2 '(R1:R2,:)'        'rmw' }
    };

  for icase = 1:Ncases
    Case = CaseList{icase};

    fprintf('*****************************************************************\n');
    fprintf('Generating storm measurements:\n');
    fprintf('  Case: %s\n', Case);
    fprintf('\n');

    % Write out measurement
    OutFile = sprintf('%s/storm_meas_tseries_%s.h5', Ddir, Case);
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
      Ndims    = MeasList{imeas}{4};
      Mselect  = MeasList{imeas}{5};
      Mmeas    = MeasList{imeas}{6};

      Mfile = regexprep(Mfile, '<CASE>', Case);

      fprintf('  Measurement: %s\n', Mname);
      fprintf('    Reading: %s (%s)\n', Mfile, Mvname);
      fprintf('    Number of dimensions: %d\n', Ndims);
      fprintf('    Selection spec: %s\n', Mselect);
      fprintf('    Measurement: %s\n', Mmeas);
      fprintf('\n');

      % MDATA will be either 1D, 2D or 3D
      %    1D: (t)
      %    2D: (r,t)
      %    3D: (r,z,t)
      MDATA = squeeze(h5read(Mfile, Mvname));
      X = squeeze(h5read(Mfile, '/x_coords'));
      Y = squeeze(h5read(Mfile, '/y_coords'));
      Z = squeeze(h5read(Mfile, '/z_coords'));
      T = squeeze(h5read(Mfile, '/t_coords'));

      Nx = length(X);
      Ny = length(Y);
      Nz = length(Z);
      Nt = length(T);

      % Apply the data selection
      SelCmd = sprintf('squeeze(MDATA%s)', Mselect);
      MDATA = eval(SelCmd);

      % Reduce data to 1D time series using Mmeas.
      %  If data is 1D, then done
      %  If data is 2D, run Mmeas on 1st dimension
      %  If data is 3D, run Mmeas on 1st, 2nd dimensions
      if (Ndims > 1)
        Tseries = zeros([ 1 Nt ]);
        for it = 1:Nt
          if (ndims(MDATA) == 2)
             % 2D -> (r,t)
             %   M_TEST <- (r)
             M_TEST = squeeze(MDATA(:,it));
          elseif (ndims(MDATA) == 3)
             % 3D -> (r,z,t)
             %   M_TEST <- (r,z)
             M_TEST = squeeze(MDATA(:,:,it));
          end

          if (strcmp(Mmeas, 'min'))
            Tseries(it) = min(M_TEST(:));
          elseif (strcmp(Mmeas, 'max') || strcmp(Mmeas, 'rmw'))
            % record the index and value of the maximum
            % if doing rmw, use the index to look up the radius
            %
            [ Val Ind ] = max(M_TEST(:));

            if (strcmp(Mmeas, 'max'))
              Tseries(it) = Val(1);
            elseif (strcmp(Mmeas, 'rmw'))
              [ Rind Zind ] = ind2sub(size(M_TEST), Ind(1));
              Tseries(it) = X(Rind);
            end
          end
        end
      else
        Tseries = MDATA;
      end

      % Write out each measurement (time series)
      Ovname = sprintf('/%s', Mname);
      h5create(OutFile, Ovname, size(Tseries));
      h5write (OutFile, Ovname, Tseries);
    end

    % Write out time values once per file
    h5create(OutFile, '/t_coords', size(T));
    h5write (OutFile, '/t_coords', T);
  end
end 
