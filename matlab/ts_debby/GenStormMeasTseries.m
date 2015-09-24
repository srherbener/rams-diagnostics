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
    { 'min_slp'  'DIAGS/hist_meas_press_<CASE>.h5' '/avg_sea_press' 0 250 0 1 'min' }

    { 'max_wind'     'DIAGS/hist_meas_speed_<CASE>.h5' '/avg_speed_t'  0 250 0 1 'max' }
    { 'max_wind_10m' 'DIAGS/hist_meas_speed_<CASE>.h5' '/avg_speed10m' 0 250 0 1 'max' }

    { 'ike'             'TsAveragedData/horiz_ke_<CASE>.h5' '/horiz_ke' 0 250 0 1 'na' }

    { 'rmw'          'DIAGS/hist_meas_speed_<CASE>.h5'  '/avg_speed_t'  0 250 0 1 'rmw' }
    { 'rmw_10m'      'DIAGS/hist_meas_speed_<CASE>.h5'  '/avg_speed10m' 0 250 0 1 'rmw' }

    { 'pcprate'      'DIAGS/hist_meas_pcprate_<CASE>.h5'  '/avg_pcprate'  0 250 0 1 'max' }
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
      Rmin     = MeasList{imeas}{4} .* 1000;  % convert km to m
      Rmax     = MeasList{imeas}{5} .* 1000;  % convert km to m
      Zmin     = MeasList{imeas}{6} .* 1000;  % convert km to m
      Zmax     = MeasList{imeas}{7} .* 1000;  % convert km to m
      Mmeas    = MeasList{imeas}{8};

      Mfile = regexprep(Mfile, '<CASE>', Case);

      fprintf('  Measurement: %s\n', Mname);
      fprintf('    Reading: %s (%s)\n', Mfile, Mvname);
      fprintf('    Rmin: %.2f\n', Rmin);
      fprintf('    Rmax: %.2f\n', Rmax);
      fprintf('    Zmin: %.2f\n', Zmin);
      fprintf('    Zmax: %.2f\n', Zmax);
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
      Z1 = find(Z >= Zmin, 1, 'first');
      Z2 = find(Z <= Zmax, 1, 'last');

      % Grab the number of dimensions. When you have a vector (size is [ 1 n ] or [ n 1 ]),
      % Ndims will be set to 2 when you really want it to be 1.
      Ndims = ndims(MDATA);
      if (Ndims == 2)
        if ((size(MDATA,1) == 1) || (size(MDATA,2) == 1))
          Ndims = 1
        end
      end

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
             M_TEST = squeeze(MDATA(R1:R2,it));
          elseif (ndims(MDATA) == 3)
             % 3D -> (r,z,t)
             %   M_TEST <- (r,z)
             M_TEST = squeeze(MDATA(R1:R2,Z1:Z2,it));
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
              Tseries(it) = R(Rind);
            end
          end
        end
      else
        Tseries = MDATA;
      end

      % Smooth out RMW time series
      if (strcmp(Mmeas, 'rmw'))
        Tseries = SmoothFillTseries(Tseries, Nt, 5);
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
