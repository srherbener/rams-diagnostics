function [ ] = GenHistOpenSpace( ConfigFile )
% GenHistOpenSpace generate histogram data for open space (between convective regions)
%
% This routine will take histogram files generated by azavg and form a new histogram
% file that represents the open space between convective regions. The input files need
% to have histogram file for all points, and a histogram file for points selected in
% convective regions. Then the output file will just be the "all points" file minus
% the "convective region selection" file.
%

% For now define mappings for what we want converted here
FprefixAll = {
  'hist_aggr_all'
  'hist_graup_all'
  'hist_hail_all'
  'hist_pris_all'
  'hist_snow_all'
  'hist_lh_tott_lht1p0'
  'hist_vaptott_vt0p5'
  };

FprefixConv = {
  'hist_aggr_twp4'
  'hist_graup_twp4'
  'hist_hail_twp4'
  'hist_pris_twp4'
  'hist_snow_twp4'
  'hist_lh_tott_twp4_lht1p0'
  'hist_vaptott_twp4_vt0p5'
  };

FprefixOpen = {
  'hist_aggr_open'
  'hist_graup_open'
  'hist_hail_open'
  'hist_pris_open'
  'hist_snow_open'
  'hist_lh_tott_open_lht1p0'
  'hist_vaptott_open_vt0p5'
  };

Vars = {
  'aggr'
  'graup'
  'hail'
  'pris'
  'snow'
  'lh_tott'
  'vaptott'
  };


% Read in the config data
[ Config ] = ReadConfig(ConfigFile);

AzavgDir = Config.AzavgDir;

for icase = 1:length(Config.Cases)
  Case = Config.Cases(icase).Cname;
  for ivar = 1: length(Vars)
    Vname = Vars{ivar};
    FpAll = FprefixAll{ivar};
    FpConv = FprefixConv{ivar};
    FpOpen = FprefixOpen{ivar};

    fprintf('***********************************************************************\n');
    fprintf('Generating open space histograms:\n');
    fprintf('  Variable: %s\n', Vname);
    fprintf('  Case: %s\n', Case);
    fprintf('\n');

    InFileAll  = sprintf('%s/%s_%s.h5', AzavgDir, FpAll, Case);
    InFileConv = sprintf('%s/%s_%s.h5', AzavgDir, FpConv, Case);
    OutFile    = sprintf('%s/%s_%s.h5', AzavgDir, FpOpen, Case);
    Hdset      = sprintf('/%s', Vname);

    % Read in the histogram data. HIST will be organized as (r,b,z,t) where
    %    x --> radius
    %    y --> histogram bins
    %    z --> heights
    %    t --> time
    fprintf('Reading file: %s\n', InFileAll);
    fprintf('\n');
    X = hdf5read(InFileAll, '/x_coords');
    Y = hdf5read(InFileAll, '/y_coords');
    Z = hdf5read(InFileAll, '/z_coords');
    T = hdf5read(InFileAll, '/t_coords');
    HIST_ALL = hdf5read(InFileAll, Hdset);

    fprintf('Reading file: %s\n', InFileConv);
    fprintf('\n');
    HIST_CONV = hdf5read(InFileConv, Hdset);

    HIST = HIST_ALL - HIST_CONV;

    % Write out the difference histogram, plus the coordinate data
    %   (make this so that GenHistMeas and GenProfMeas can use it.)
    fprintf('Writing file: %s\n', OutFile);
    fprintf('\n');

    hdf5write(OutFile, Vname, HIST);

    hdf5write(OutFile, '/x_coords', X, 'WriteMode', 'append');
    hdf5write(OutFile, '/y_coords', Y, 'WriteMode', 'append');
    hdf5write(OutFile, '/z_coords', Z, 'WriteMode', 'append');
    hdf5write(OutFile, '/t_coords', T, 'WriteMode', 'append');

    fprintf('\n');
  end
end

end
