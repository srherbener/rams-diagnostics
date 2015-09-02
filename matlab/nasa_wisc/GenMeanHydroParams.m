function [ ] = GenMeanHydroParams()
% GenMeanHydroParams generate mean hydrometeor parameters for one-moment sim
%                    based on data from two-moment sim

  Ddir = 'DIAGS';

  % make sure output directory exists
  if (exist(Ddir, 'dir') ~= 7)
    mkdir(Ddir);
  end

%  Case = 'RCE_S300';  % 2 moment case
  Case = 'RCE_S300_UB5';  % 2 moment case with stronger surface interaction

  % input_file input_dataset output_file output_dataset hydrometeor_index apply_volume_density_weighting
  VarSets = ...
    {
      % CLOUD
      { 'HDF5/TsAveragedData/hda_cloud_0p01_<CASE>.h5'      '/hda_cloud'      'DIAGS/avg_cloud_<CASE>.h5'      '/avg_cloud'      1 'hda' }
      { 'HDF5/TsAveragedData/hda_cloud_num_0p01_<CASE>.h5'  '/hda_cloud_num'  'DIAGS/avg_cloud_num_<CASE>.h5'  '/avg_cloud_num'  1 'hda' }
      { 'HDF5/TsAveragedData/hda_cloud_diam_0p01_<CASE>.h5' '/hda_cloud_diam' 'DIAGS/avg_cloud_diam_<CASE>.h5' '/avg_cloud_diam' 0 'hda' }

      { 'HDF5/TsAveragedData/hist_cloud_0p01_<CASE>.h5'      '/hist_cloud'      'DIAGS/havg_cloud_<CASE>.h5'      '/avg_cloud'      1 'hist_log' }
      { 'HDF5/TsAveragedData/hist_cloud_num_0p01_<CASE>.h5'  '/hist_cloud_num'  'DIAGS/havg_cloud_num_<CASE>.h5'  '/avg_cloud_num'  1 'hist_log' }
      { 'HDF5/TsAveragedData/hist_cloud_diam_0p01_<CASE>.h5' '/hist_cloud_diam' 'DIAGS/havg_cloud_diam_<CASE>.h5' '/avg_cloud_diam' 0 'hist_log' }
  
      % DRIZZLE
      { 'HDF5/TsAveragedData/hda_driz_0p01_<CASE>.h5'      '/hda_driz'      'DIAGS/avg_driz_<CASE>.h5'      '/avg_driz'      1 'hda' }
      { 'HDF5/TsAveragedData/hda_driz_num_0p01_<CASE>.h5'  '/hda_driz_num'  'DIAGS/avg_driz_num_<CASE>.h5'  '/avg_driz_num'  1 'hda' }
      { 'HDF5/TsAveragedData/hda_driz_diam_0p01_<CASE>.h5' '/hda_driz_diam' 'DIAGS/avg_driz_diam_<CASE>.h5' '/avg_driz_diam' 0 'hda' }
  
      { 'HDF5/TsAveragedData/hist_driz_0p01_<CASE>.h5'      '/hist_driz'      'DIAGS/havg_driz_<CASE>.h5'      '/avg_driz'      1 'hist_log' }
      { 'HDF5/TsAveragedData/hist_driz_num_0p01_<CASE>.h5'  '/hist_driz_num'  'DIAGS/havg_driz_num_<CASE>.h5'  '/avg_driz_num'  1 'hist_log' }
      { 'HDF5/TsAveragedData/hist_driz_diam_0p01_<CASE>.h5' '/hist_driz_diam' 'DIAGS/havg_driz_diam_<CASE>.h5' '/avg_driz_diam' 0 'hist_log' }
  
      % RAIN
      { 'HDF5/TsAveragedData/hda_rain_0p01_<CASE>.h5'      '/hda_rain'      'DIAGS/avg_rain_<CASE>.h5'      '/avg_rain'      1 'hda' }
      { 'HDF5/TsAveragedData/hda_rain_num_0p01_<CASE>.h5'  '/hda_rain_num'  'DIAGS/avg_rain_num_<CASE>.h5'  '/avg_rain_num'  1 'hda' }
      { 'HDF5/TsAveragedData/hda_rain_diam_0p01_<CASE>.h5' '/hda_rain_diam' 'DIAGS/avg_rain_diam_<CASE>.h5' '/avg_rain_diam' 0 'hda' }
  
      { 'HDF5/TsAveragedData/hist_rain_0p01_<CASE>.h5'      '/hist_rain'      'DIAGS/havg_rain_<CASE>.h5'      '/avg_rain'      1 'hist_log' }
      { 'HDF5/TsAveragedData/hist_rain_num_0p01_<CASE>.h5'  '/hist_rain_num'  'DIAGS/havg_rain_num_<CASE>.h5'  '/avg_rain_num'  1 'hist_log' }
      { 'HDF5/TsAveragedData/hist_rain_diam_0p01_<CASE>.h5' '/hist_rain_diam' 'DIAGS/havg_rain_diam_<CASE>.h5' '/avg_rain_diam' 0 'hist_log' }
  
      % PRIS
      { 'HDF5/TsAveragedData/hda_pris_0p01_<CASE>.h5'      '/hda_pris'      'DIAGS/avg_pris_<CASE>.h5'      '/avg_pris'      1 'hda' }
      { 'HDF5/TsAveragedData/hda_pris_num_0p01_<CASE>.h5'  '/hda_pris_num'  'DIAGS/avg_pris_num_<CASE>.h5'  '/avg_pris_num'  1 'hda' }
      { 'HDF5/TsAveragedData/hda_pris_diam_0p01_<CASE>.h5' '/hda_pris_diam' 'DIAGS/avg_pris_diam_<CASE>.h5' '/avg_pris_diam' 0 'hda' }
  
      { 'HDF5/TsAveragedData/hist_pris_0p01_<CASE>.h5'      '/hist_pris'      'DIAGS/havg_pris_<CASE>.h5'      '/avg_pris'      1 'hist_log' }
      { 'HDF5/TsAveragedData/hist_pris_num_0p01_<CASE>.h5'  '/hist_pris_num'  'DIAGS/havg_pris_num_<CASE>.h5'  '/avg_pris_num'  1 'hist_log' }
      { 'HDF5/TsAveragedData/hist_pris_diam_0p01_<CASE>.h5' '/hist_pris_diam' 'DIAGS/havg_pris_diam_<CASE>.h5' '/avg_pris_diam' 0 'hist_log' }
  
      % SNOW
      { 'HDF5/TsAveragedData/hda_snow_0p01_<CASE>.h5'      '/hda_snow'      'DIAGS/avg_snow_<CASE>.h5'      '/avg_snow'      1 'hda' }
      { 'HDF5/TsAveragedData/hda_snow_num_0p01_<CASE>.h5'  '/hda_snow_num'  'DIAGS/avg_snow_num_<CASE>.h5'  '/avg_snow_num'  1 'hda' }
      { 'HDF5/TsAveragedData/hda_snow_diam_0p01_<CASE>.h5' '/hda_snow_diam' 'DIAGS/avg_snow_diam_<CASE>.h5' '/avg_snow_diam' 0 'hda' }
  
      { 'HDF5/TsAveragedData/hist_snow_0p01_<CASE>.h5'      '/hist_snow'      'DIAGS/havg_snow_<CASE>.h5'      '/avg_snow'      1 'hist_log' }
      { 'HDF5/TsAveragedData/hist_snow_num_0p01_<CASE>.h5'  '/hist_snow_num'  'DIAGS/havg_snow_num_<CASE>.h5'  '/avg_snow_num'  1 'hist_log' }
      { 'HDF5/TsAveragedData/hist_snow_diam_0p01_<CASE>.h5' '/hist_snow_diam' 'DIAGS/havg_snow_diam_<CASE>.h5' '/avg_snow_diam' 0 'hist_log' }
  
      % AGGREGATES
      { 'HDF5/TsAveragedData/hda_aggr_0p01_<CASE>.h5'      '/hda_aggr'      'DIAGS/avg_aggr_<CASE>.h5'      '/avg_aggr'      1 'hda' }
      { 'HDF5/TsAveragedData/hda_aggr_num_0p01_<CASE>.h5'  '/hda_aggr_num'  'DIAGS/avg_aggr_num_<CASE>.h5'  '/avg_aggr_num'  1 'hda' }
      { 'HDF5/TsAveragedData/hda_aggr_diam_0p01_<CASE>.h5' '/hda_aggr_diam' 'DIAGS/avg_aggr_diam_<CASE>.h5' '/avg_aggr_diam' 0 'hda' }
  
      { 'HDF5/TsAveragedData/hist_aggr_0p01_<CASE>.h5'      '/hist_aggr'      'DIAGS/havg_aggr_<CASE>.h5'      '/avg_aggr'      1 'hist_log' }
      { 'HDF5/TsAveragedData/hist_aggr_num_0p01_<CASE>.h5'  '/hist_aggr_num'  'DIAGS/havg_aggr_num_<CASE>.h5'  '/avg_aggr_num'  1 'hist_log' }
      { 'HDF5/TsAveragedData/hist_aggr_diam_0p01_<CASE>.h5' '/hist_aggr_diam' 'DIAGS/havg_aggr_diam_<CASE>.h5' '/avg_aggr_diam' 0 'hist_log' }
  
      % GRAUPEL
      { 'HDF5/TsAveragedData/hda_graup_0p01_<CASE>.h5'      '/hda_graup'      'DIAGS/avg_graup_<CASE>.h5'      '/avg_graup'      1 'hda' }
      { 'HDF5/TsAveragedData/hda_graup_num_0p01_<CASE>.h5'  '/hda_graup_num'  'DIAGS/avg_graup_num_<CASE>.h5'  '/avg_graup_num'  1 'hda' }
      { 'HDF5/TsAveragedData/hda_graup_diam_0p01_<CASE>.h5' '/hda_graup_diam' 'DIAGS/avg_graup_diam_<CASE>.h5' '/avg_graup_diam' 0 'hda' }
  
      { 'HDF5/TsAveragedData/hist_graup_0p01_<CASE>.h5'      '/hist_graup'      'DIAGS/havg_graup_<CASE>.h5'      '/avg_graup'      1 'hist_log' }
      { 'HDF5/TsAveragedData/hist_graup_num_0p01_<CASE>.h5'  '/hist_graup_num'  'DIAGS/havg_graup_num_<CASE>.h5'  '/avg_graup_num'  1 'hist_log' }
      { 'HDF5/TsAveragedData/hist_graup_diam_0p01_<CASE>.h5' '/hist_graup_diam' 'DIAGS/havg_graup_diam_<CASE>.h5' '/avg_graup_diam' 0 'hist_log' }
  
      % HAIL
      { 'HDF5/TsAveragedData/hda_hail_0p01_<CASE>.h5'      '/hda_hail'      'DIAGS/avg_hail_<CASE>.h5'      '/avg_hail'      1 'hda' }
      { 'HDF5/TsAveragedData/hda_hail_num_0p01_<CASE>.h5'  '/hda_hail_num'  'DIAGS/avg_hail_num_<CASE>.h5'  '/avg_hail_num'  1 'hda' }
      { 'HDF5/TsAveragedData/hda_hail_diam_0p01_<CASE>.h5' '/hda_hail_diam' 'DIAGS/avg_hail_diam_<CASE>.h5' '/avg_hail_diam' 0 'hda' }

      { 'HDF5/TsAveragedData/hist_hail_0p01_<CASE>.h5'      '/hist_hail'      'DIAGS/havg_hail_<CASE>.h5'      '/avg_hail'      1 'hist_log' }
      { 'HDF5/TsAveragedData/hist_hail_num_0p01_<CASE>.h5'  '/hist_hail_num'  'DIAGS/havg_hail_num_<CASE>.h5'  '/avg_hail_num'  1 'hist_log' }
      { 'HDF5/TsAveragedData/hist_hail_diam_0p01_<CASE>.h5' '/hist_hail_diam' 'DIAGS/havg_hail_diam_<CASE>.h5' '/avg_hail_diam' 0 'hist_log' }
  
    };
  Nset = length(VarSets);

  % reference density taken from RAMS output
  Dens = [
    1.152
    1.146
    1.140
    1.133
    1.126
    1.119
    1.111
    1.102
    1.093
    1.083
    1.072
    1.062
    1.050
    1.038
    1.026
    1.012
    0.998
    0.983
    0.969
    0.953
    0.937
    0.920
    0.903
    0.885
    0.867
    0.847
    0.828
    0.809
    0.788
    0.767
    0.745
    0.724
    0.700
    0.675
    0.651
    0.628
    0.603
    0.578
    0.553
    0.528
    0.503
    0.478
    0.453
    0.429
    0.406
    0.385
    0.364
    0.344
    0.326
    0.307
    0.290
    0.273
    0.257
    0.241
    0.226
    0.212
    0.198
    0.183
    0.164
    0.148
    0.133
    0.122
    0.112
    0.102
    0.094
    0.086
    0.079
    0.072
    0.066
    0.060
    0.055
    0.050
    0.046
    0.042
    0.039
    ];

  fprintf('***************************************************************\n');
  fprintf('Generating parameter averages:\n');
  fprintf('  Case: %s\n', Case);
  fprintf('\n');


  for iset = 1:Nset
    InFile     = regexprep(VarSets{iset}{1}, '<CASE>', Case);
    InVname    = VarSets{iset}{2};
    OutFile    = regexprep(VarSets{iset}{3}, '<CASE>', Case);
    OutVname   = VarSets{iset}{4};
    UseWeights = VarSets{iset}{5};
    InType     = VarSets{iset}{6};

    % Read in averages, format is (2,z,t) where
    %   (1,z,t) are the sums (per level)
    %   (2,z,t) are the counts (per level)
    fprintf('  Reading (%s): %s (%s)\n', InType, InFile, InVname);
    fprintf('\n');

    HDATA = squeeze(h5read(InFile, InVname));
    X     = squeeze(h5read(InFile, '/x_coords'));
    Y     = squeeze(h5read(InFile, '/y_coords'));
    Z     = squeeze(h5read(InFile, '/z_coords'));
    T     = squeeze(h5read(InFile, '/t_coords'));
    Nx    = length(X);
    Nz    = length(Z);
    Nt    = length(T);

    % Select all time points for now.
    T1 = 1;
%    T1 = 50;
    T2 = Nt;

    % Get stretched vertical grid spacing
    DeltaZ = Z(2:end) - Z(1:end-1);
    DeltaZ = [ DeltaZ(1) DeltaZ' ]';  % just repeat the lower delta for the first entry

    % Mixing ratio is kg/kg(air) and number concentration is #/kg(air). The horizontal grid
    % space is constant everywhere, but the vertical grid space varies. Also the density
    % varies with height. In order to account for this create weights for averaging that
    % are based on the relative amounts of air at each vertical level. The volume changes by
    % a factor of DeltaZ(k)/DeltaZ(1) at each level, and the density changes by a factor of
    % Dens(k)/Dens(1) at each level. The overall factor is then:
    %
    %     Weights(k) = [ DeltaZ(k) * Dens(k) ] / [ DeltaZ(1) * Dens(1) ]
    %
    if (UseWeights == 1)
      % using weights
      WEIGHTS = (DeltaZ .* Dens) ./ (DeltaZ(1) * Dens(1));
    else
      % not using weights, set WEIGHTS to all ones
      WEIGHTS = ones([ Nz 1 ]);
    end

    % Calculate the averages. Understand three types of input of which all come
    % runnign tsavg.
    %
    %   'hda' -> horizontal domain average
    %   'hist_lin' -> histogram with linear bin edges
    %   'hist_log' -> histogram with logrithmic bin edges
    %
    if (strcmp(InType, 'hda'))
      % Horzontal domain averages

      % Input data is organized as (2,z,t), where first dimension is:
      %     (1,z,t) -> sums
      %     (2,z,t) -> counts

      % extract sum and count values, also select time range for averaging
      SUM   = squeeze(HDATA(1,:,T1:T2));    % sum 
      COUNT = squeeze(HDATA(2,:,T1:T2));  % count

      % SUM and COUNT are both (z,t), sum up both of these across
      % the time dimension yielding vectors of length Nz.
      SUM   = sum(SUM,2);
      COUNT = sum(COUNT,2);
  
      % The volume of grid cells and the density of air are constant for a given
      % level, so don't apply the weights to get the profiles. But do apply weights
      % for total domain averages.
      %
      % For total domain averages: divide the weighted average of SUM across all levels
      % by weighted average of COUNT across all levels. This yields:
      %
      %    { sum(SUM * WEIGHTS) / sum( WEIGHTS) } / { sum(COUNT * WEIGHTS) / sum(WEIGHTS) }
      %
      % The sum(WEIGHTS) terms cancel, leaving:
      %
      %     sum(SUM * WEIGHTS) / sum(COUNT * WEIGHTS)
      %
      AVG_PROF = SUM ./ COUNT;
      AVG_DOM = squeeze(sum(SUM .* WEIGHTS)) / squeeze(sum(COUNT .* WEIGHTS));
    elseif (strncmp(InType, 'hist', 4))
      % Histograms

      % Input data is organized as (x,z,t) where x contains the bin values
      %
      % If have logrithmic bin values, do the averaging on the logrithm of the
      % bin values and after the averaging take the anit-logrithm for the result.
      %
      if (strcmp(InType, 'hist_log'))
        BINS = log10(X);
      elseif (strcmp(InType, 'hist_lin'))
        BINS = X;
      end

      % Sum up counts along the time dimension to get a vertical profile of histograms.
      % Then do a weighted average along the vertical to get a single domain histogram.
      %
      % Use the fractional area with a parameter of 0.5 when reducing the histograms.
      % This will produce the median value.
      HIST_WEIGHTS = repmat(WEIGHTS', [ Nx 1 ]);
      HIST_PROF = sum(HDATA,3);
      HIST = sum((HIST_PROF .* HIST_WEIGHTS), 2) ./ sum(HIST_WEIGHTS, 2);

      AVG_PROF = ReduceHists(HIST_PROF, 1, BINS, 'farea', 0.5);
      AVG_DOM  = ReduceHists(HIST, 1, BINS, 'farea', 0.5);

      % If the input type is 'hist_log', then the results in AVG_PROF and AVG_DOM
      % are base 10 logarithm values. Convert these back to actual values.
      if (strcmp(InType, 'hist_log'))
        AVG_PROF = 10 .^ AVG_PROF;
        AVG_DOM  = 10 ^ AVG_DOM;
      end
      
    end

    % Output
    OutProfVname = sprintf('%s_prof', OutVname);

    fprintf('  Writing: %s (%s)\n', OutFile, OutVname);
    fprintf('  Writing: %s (%s)\n', OutFile, OutProfVname);
    fprintf('\n');

    % Remove existing file so that subsequent create/write command will not
    % have to deal with replacing data.
    if (exist(OutFile, 'file') == 2)
      delete(OutFile);
    end

    h5create(OutFile, OutVname,     size(AVG_DOM));
    h5create(OutFile, OutProfVname, size(AVG_PROF));
    h5create(OutFile, '/z_coords',  size(Z));

    h5write(OutFile, OutVname,     AVG_DOM);
    h5write(OutFile, OutProfVname, AVG_PROF);
    h5write(OutFile, '/z_coords',  Z);
  end
end
