function [ ] = GenTopSwup(ConfigFile)
% GenTopSwup generate time series of the top level average SW upward flux

% Read the config file to get the structure of how the data is laid out in
% the file system.
[ Config ] = ReadConfig(ConfigFile);

% This is a simple script that assumes that you already have the domain averge CWP and RWP
% in hdf5 files.
    
Tdir = Config.TsavgDir;
SwupVar = 'hda_swup';
OutVar  = 'hda_top_swup';

% Find and replace underscores in Ptitle, Ylabel with blank spaces
for icase = 1:length(Config.Cases)
    Case = Config.Cases(icase).Cname;
    SwupFile = sprintf('%s/%s_%s.h5', Tdir, SwupVar, Case);
    OutFile = sprintf('%s/%s_%s.h5', Tdir, OutVar, Case);

    fprintf('***************************************************************\n');
    fprintf('Extracting top level SW upward flux:\n');
    fprintf('  Case: %s\n', Case);
    fprintf('  Input swup file: %s\n', SwupFile);
    fprintf('  Output swup file: %s\n', OutFile);
    fprintf('\n');

    % read in swup, (x,y,z,t)
    SWUP  = hdf5read(SwupFile, SwupVar);

    % copy the coordinate values to the output to keep GenTsPlots happy
    Xcoords = hdf5read(SwupFile, 'x_coords');
    Ycoords = hdf5read(SwupFile, 'y_coords');
    Zcoords = hdf5read(SwupFile, 'z_coords');
    Tcoords = hdf5read(SwupFile, 't_coords');

    % RAMS/REVU set the top boundary to all zeros so take the top-1 level

    SWUP_TOP = SWUP(:,:,end-1,:);
    ZcoordsTop = Zcoords(end-1);

    hdf5write(OutFile, OutVar, SWUP_TOP);
    hdf5write(OutFile, 'x_coords', Xcoords, 'WriteMode', 'append');
    hdf5write(OutFile, 'y_coords', Ycoords, 'WriteMode', 'append');
    hdf5write(OutFile, 'z_coords', ZcoordsTop, 'WriteMode', 'append');
    hdf5write(OutFile, 't_coords', Tcoords, 'WriteMode', 'append');

end
