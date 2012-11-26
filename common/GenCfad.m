function [ ] = GenCfad( ConfigFile )
% GenCfad generate CFAD data
%
% This routine will select histogram data from given radius range
% and time range. Then the histogram count data will simply be
% summed up over the selected space. This will result in 2D
% data of the form: (b,z) where b is the histogram bins and
% z is height.

% Read in the config data
[ Config ] = ReadConfig(ConfigFile);

DiagDir = Config.DiagDir;

% Make sure output directory exists
if (exist(DiagDir, 'dir') ~= 7)
    mkdir(DiagDir);
end

for icase = 1:length(Config.Cases)
  Case = Config.Cases(icase).Cname;
  for icfad = 1: length(Config.Cfad)
    Name = Config.Cfad(icfad).Name;
    InDir = Config.Cfad(icfad).InDir;
    Fprefix = Config.Cfad(icfad).Fprefix;
    Vname = Config.Cfad(icfad).Rvar;
    Rmin = Config.Cfad(icfad).Rmin;
    Rmax = Config.Cfad(icfad).Rmax;
    Tmin = Config.Cfad(icfad).Tmin;
    Tmax = Config.Cfad(icfad).Tmax;

    fprintf('***********************************************************************\n');
    fprintf('Generating CFAD:\n');
    fprintf('  Name: %s\n', Name);
    fprintf('  Variable: %s\n', Vname);
    fprintf('  Case: %s\n', Case);
    fprintf('  Rmin, Rmax: %.2f, %.2f\n', Rmin, Rmax);
    fprintf('  Tmin, Tmax: %.2f, %.2f\n', Tmin, Tmax);
    fprintf('\n');

    InFile  = sprintf('%s/%s_%s.h5', InDir, Fprefix, Case);
    OutFile = sprintf('%s/%s_%s.h5', DiagDir, Name, Case);
    Hdset   = sprintf('/%s', Vname);

    % Read in the histogram data. HIST will be organized as (r,b,z,t) where
    %    r --> radius
    %    b --> histogram bins
    %    z --> heights
    %    t --> time
    fprintf('Reading file: %s\n', InFile);
    fprintf('\n');
    R = hdf5read(InFile, '/x_coords') / 1000; % km
    B = hdf5read(InFile, '/y_coords');
    Z = hdf5read(InFile, '/z_coords');
    T = hdf5read(InFile, '/t_coords')/ 3600; % hr
    HIST = hdf5read(InFile, Hdset);

    % Select the space to do the summing
    R1 = find(R >= Rmin, 1, 'first');
    R2 = find(R <= Rmax, 1, 'last');
    T1 = find(T >= Tmin, 1, 'first');
    T2 = find(T <= Tmax, 1, 'last');

    HSEL = HIST(R1:R2,:,:,T1:T2);

    % Sum across the r and t dimensions (1st and 4th)
    CFAD = squeeze(sum(sum(HSEL,1),4));

    fprintf('Writing file: %s\n', OutFile);
    fprintf('\n');

    Hdset = sprintf('/%s', Vname);
    hdf5write(OutFile, Hdset, CFAD);

    hdf5write(OutFile, '/Bins', B, 'WriteMode', 'append');
    hdf5write(OutFile, '/Height', Z, 'WriteMode', 'append');

    fprintf('\n');
  end
end

end
