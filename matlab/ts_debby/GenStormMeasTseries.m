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

  % Description of measurements
  MeasList = {
    { 'min_slp_fa'  'DIAGS/hist_meas_press' '/avg_sea_press_fa' '(1:28,:)'      'min' } % selection is R = 0 to 250 km
    { 'min_slp_wm'  'DIAGS/hist_meas_press' '/avg_sea_press_wm' '(1:28,:)'      'min' } % selection is R = 0 to 250 km

    { 'max_wind_fa'     'DIAGS/hist_meas_speed' '/avg_speed_t_fa'   '(1:28,2:14,:)' 'max' } % selection is R = 0 to 250 km, Z = 0 to 1000 m
    { 'max_wind_wm'     'DIAGS/hist_meas_speed' '/avg_speed_t_wm'   '(1:28,2:14,:)' 'max' } % selection is R = 0 to 250 km, Z = 0 to 1000 m
    { 'max_wind_10m_fa' 'DIAGS/hist_meas_speed' '/avg_speed10m_fa'  '(1:28,:)'      'max' } % selection is R = 0 to 250 km
    { 'max_wind_10m_wm' 'DIAGS/hist_meas_speed' '/avg_speed10m_wm'  '(1:28,:)'      'max' } % selection is R = 0 to 250 km
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
      Mfprefix = MeasList{imeas}{2};
      Mvname   = MeasList{imeas}{3};
      Mselect  = MeasList{imeas}{4};
      Mmeas    = MeasList{imeas}{5};

      Mfile = sprintf('%s_%s.h5', Mfprefix, Case);

      fprintf('  Measurement: %s\n', Mname);
      fprintf('    Reading: %s (%s)\n', Mfile, Mvname);
      fprintf('    Selection spec: %s\n', Mselect);
      fprintf('    Measurement: %s\n', Mmeas);
      fprintf('\n');

      % Read in data, which is coming from either azavg or tsavg, meaning
      % that it will be 4D -> (x,y,z,t) regardless if all four dimensions
      % are actually used. Mname will indicate what size MDATA is and what
      % each of X, Y, Z, T represent.
      %
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
      for idim = 1:ndims(MDATA)-1
        % Keep reducing 1st dimension until finished
        if (strcmp(Mmeas, 'min'))
          MDATA = squeeze(min(MDATA, [] , 1));
        elseif (strcmp(Mmeas, 'max'))
          MDATA = squeeze(max(MDATA, [] , 1));
        end
      end

      % Write out each measurement (time series)
      Ovname = sprintf('/%s', Mname);
      h5create(OutFile, Ovname, size(MDATA));
      h5write (OutFile, Ovname, MDATA);
    end

    % Write out time values once per file
    h5create(OutFile, '/t_coords', size(T));
    h5write (OutFile, '/t_coords', T);
  end
end 
