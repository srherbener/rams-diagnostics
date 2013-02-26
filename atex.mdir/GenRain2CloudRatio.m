function [ ] = GenRain2CloudRatio(ConfigFile)
% GenRain2CloudRatio generate time series of the ratio of domain average CWP to domain average RWP

% Read the config file to get the structure of how the data is laid out in
% the file system.
[ Config ] = ReadConfig(ConfigFile);

% This is a simple script that assumes that you already have the domain averge CWP and RWP
% in hdf5 files.
    
Tdir = Config.TsavgDir;
RainVar = 'hda_vint_rain';
CloudVar = 'hda_vint_cloud';
OutVar = 'hda_vint_r2c';

% Find and replace underscores in Ptitle, Ylabel with blank spaces
for icase = 1:length(Config.Cases)
    Case = Config.Cases(icase).Cname;
    RainFile = sprintf('%s/%s_%s.h5', Tdir, RainVar, Case);
    CloudFile = sprintf('%s/%s_%s.h5', Tdir, CloudVar, Case);
    OutFile = sprintf('%s/%s_%s.h5', Tdir, OutVar, Case);

    fprintf('***************************************************************\n');
    fprintf('Generating rain water path to cloud water path ratio:\n');
    fprintf('  Case: %s\n', Case);
    fprintf('  Input rain file: %s\n', RainFile);
    fprintf('  Input cloud file: %s\n', CloudFile);
    fprintf('  Output file: %s\n', OutFile);
    fprintf('\n');

    % Read in the rain and cloud data, and create the ratio: rain / cloud.
    % If cloud is a zero:
    %   If rain is a zero, set the ratio to zero
    %   If rain is non-zero, set the ratio to nan.
    %
    % Rain and cloud data are 4D, (x,y,z,t), but the x,y and z dimensions
    % are all size 1. Squeeze the dimensions out so that rain and cloud are
    % just one dimension.

    RAIN  = squeeze(hdf5read(RainFile, RainVar));
    CLOUD = squeeze(hdf5read(CloudFile, CloudVar));

    % copy the coordinate values to the output to keep GenTsPlots happy
    Xcoords = hdf5read(CloudFile, 'x_coords');
    Ycoords = hdf5read(CloudFile, 'y_coords');
    Zcoords = hdf5read(CloudFile, 'z_coords');
    Tcoords = hdf5read(CloudFile, 't_coords');

    for i = 1:length(RAIN)
      if (CLOUD(i) == 0)
        % When 0 / 0 --> allow this to be zero,
        % otherwise if n / 0 set to nan.
        if (RAIN(i) == 0)
           R2C(i) = 0;
        else
           R2C(i) = nan;
        end
      else 
        % okay to divide, but avoid taking ratio of tiny numbers (which
        % would be misleading)
        if ((abs(RAIN(i)) < 1.0e-6) && (abs(CLOUD(i) < 1.0e-6)))
          R2C(i) = 0;
        else
          R2C(i) = RAIN(i) / CLOUD(i);
        end
      end
    end

    hdf5write(OutFile, OutVar, R2C);
    hdf5write(OutFile, 'x_coords', Xcoords, 'WriteMode', 'append');
    hdf5write(OutFile, 'y_coords', Ycoords, 'WriteMode', 'append');
    hdf5write(OutFile, 'z_coords', Zcoords, 'WriteMode', 'append');
    hdf5write(OutFile, 't_coords', Tcoords, 'WriteMode', 'append');

end
