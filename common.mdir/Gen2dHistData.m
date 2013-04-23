function [ ] = Gen2dHistData(ConfigFile)
% Gen2dHistData form 2d Histogram data via data selection from tsavg diagnostics

% Read the config file to get the structure of how the data is laid out in
% the file system.
[ Config ] = ReadConfig(ConfigFile);
    
UndefVal = Config.UndefVal;

Ddir = Config.DiagDir;
% make sure output directory exists
if (exist(Ddir, 'dir') ~= 7)
  mkdir(Ddir);
end

Fsize = 45;

for icase = 1:length(Config.Cases)
  Case  = Config.Cases(icase).Cname;
  Pname = Config.Cases(icase).Pname;
  for ihist = 1:length(Config.Hist2d)
    InFile  = sprintf('%s_%s.h5', Config.Hist2d(ihist).Fprefix, Case);
    Var     = Config.Hist2d(ihist).Var;
    Xmin    = Config.Hist2d(ihist).Xmin;
    Xmax    = Config.Hist2d(ihist).Xmax;
    Xgroup  = Config.Hist2d(ihist).Xgroup;
    Ymin    = Config.Hist2d(ihist).Ymin;
    Ymax    = Config.Hist2d(ihist).Ymax;
    Ygroup  = Config.Hist2d(ihist).Ygroup;
    Tmin    = Config.Hist2d(ihist).Tmin;
    Tmax    = Config.Hist2d(ihist).Tmax;
    OutFile = sprintf('%s/%s_%s.h5', Ddir, Config.Hist2d(ihist).Name, Case);

    fprintf('********************************************************************\n');
    fprintf('Generating 2D histogram data:\n');
    fprintf('  Input file: %s\n', InFile);
    fprintf('  Input variable: %s\n', Var);
    fprintf('  Data selection:\n');
    fprintf('    Xrange, Xgroup: [ %.2f %.2f], %d\n', Xmin, Xmax, Xgroup);
    fprintf('    Yrange, Ygroup: [ %.2f %.2f], %d\n', Ymin, Ymax, Ygroup);
    fprintf('    Trange: [ %.2f %.2f]\n', Tmin, Tmax);
    fprintf('  Output File: %s\n', OutFile);
    fprintf('\n');

    % Read in the counts and do the data selection (along with combining bins)
    C = hdf5read(InFile, Var);
    X = hdf5read(InFile, 'x_coords');
    Y = hdf5read(InFile, 'y_coords');
    Z = hdf5read(InFile, 'z_coords');
    T = hdf5read(InFile, 't_coords')/3600; % hr
 
    T1 = find(T >= Tmin, 1, 'first');
    T2 = find(T <= Tmax, 1, 'last');
 
    [ COUNTS, XL, XU, YL, YU ] = GenCountBins(C(:,:,:,T1:T2), X, Y, Xmin, Xmax, Xgroup, Ymin, Ymax, Ygroup);
    TIMES = T(T1:T2);

    % If the time range was more than a single point, then sum up counts across the time
    % dimension. Time will be the last dimension.
    if ((T2-T1) > 0)
      COUNTS = sum(COUNTS, ndims(COUNTS));
    end

    % output
    fprintf('Writing file: %s\n', OutFile);
    fprintf('\n');

    hdf5write(OutFile, '/COUNTS', COUNTS);

    hdf5write(OutFile, '/XL',    XL, 'WriteMode', 'append');
    hdf5write(OutFile, '/XU',    XU, 'WriteMode', 'append');
    hdf5write(OutFile, '/YL',    YL, 'WriteMode', 'append');
    hdf5write(OutFile, '/YU',    YU, 'WriteMode', 'append');
    hdf5write(OutFile,  '/T', TIMES, 'WriteMode', 'append');

    fprintf('\n');
  end
end

end
