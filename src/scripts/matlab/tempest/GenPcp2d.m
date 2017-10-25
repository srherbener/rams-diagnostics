function [] = GenPcp2d(ConfigFile)
% GenPcp2d generate time series of accum precip and precip rate for 2D sims

    % Read the config file to get the structure of how the data is laid out in
    % the file system.
    [ Config ] = ReadConfig(ConfigFile);

    Ddir = Config.DiagDir;
    Hdir = 'HDF5';

    if (exist(Ddir, 'dir') ~= 7)
      mkdir(Ddir);
    end

    AccumPrecipName = 'accpr';
    AccumPrecipFprefix = 'accpr';
    PrecipRateName = 'pcprr';
    PrecipRateFprefix = 'pcprr';

    for icase = 1:length(Config.Cases)
        Case = Config.Cases(icase).Cname;
        AccumPrecipFile = sprintf('%s/%s/%s-a-AS-1999-02-10-040000-g1.h5', Hdir, Case, AccumPrecipFprefix);
        PrecipRateFile  = sprintf('%s/%s/%s-a-AS-1999-02-10-040000-g1.h5', Hdir, Case, PrecipRateFprefix);

        fprintf('***************************************************************\n');
        fprintf('Generating 2D precip data:\n');
        fprintf('  Case: %s\n', Case);
        fprintf('  Accum precip file: %s\n', AccumPrecipFile);
        fprintf('    Var name: %s\n', AccumPrecipName);
        fprintf('  Precip rate file: %s\n', PrecipRateFile);
        fprintf('    Var name: %s\n', PrecipRateName);
        fprintf('\n');
 
        % data will be organized as (x,t)
        %   x is zonal dimension
        %   t is time
        ACCPR = squeeze(hdf5read(AccumPrecipFile, AccumPrecipName));
        PCPRR = squeeze(hdf5read(PrecipRateFile, PrecipRateName));
        T   = squeeze(hdf5read(AccumPrecipFile, 't_coords'));

        % For accum precip, just som over the domain
        % For precip rate, create mean over domain with and without zeros
        ACC_PCP = squeeze(sum(ACCPR,1));

        PCP_RATE = squeeze(mean(PCPRR,1));
        PCPRR_NZ = PCPRR;
        PCPRR_NZ(PCPRR_NZ == 0) = nan;
        PCP_RATE_NZ = squeeze(nanmean(PCPRR_NZ,1));

        % Output coordinates
        X = 1;
        Y = 1;
        Z = 1;

        Nx = 1;
        Ny = 1;
        Nz = 1;
        Nt = length(T);

        % write out accum precip
        OutFile = sprintf('%s/accpr_%s.h5', Ddir, Case);
        fprintf('Writing: %s\n', OutFile);

        Ovar = reshape(ACC_PCP, [ Nx Ny Nz Nt ]);

        hdf5write(OutFile, '/accpr', Ovar);
        hdf5write(OutFile, 'x_coords', X, 'WriteMode', 'append');
        hdf5write(OutFile, 'y_coords', Y, 'WriteMode', 'append');
        hdf5write(OutFile, 'z_coords', Z, 'WriteMode', 'append');
        hdf5write(OutFile, 't_coords', T, 'WriteMode', 'append');
        fprintf('\n');

        % write out precip rate
        OutFile = sprintf('%s/pcprr_%s.h5', Ddir, Case);
        fprintf('Writing: %s\n', OutFile);

        Ovar = reshape(PCP_RATE, [ Nx Ny Nz Nt ]);
        hdf5write(OutFile, '/pcprr', Ovar);

        Ovar = reshape(PCP_RATE_NZ, [ Nx Ny Nz Nt ]);
        hdf5write(OutFile, '/pcprr_nz', Ovar, 'WriteMode', 'append');

        hdf5write(OutFile, 'x_coords', X, 'WriteMode', 'append');
        hdf5write(OutFile, 'y_coords', Y, 'WriteMode', 'append');
        hdf5write(OutFile, 'z_coords', Z, 'WriteMode', 'append');
        hdf5write(OutFile, 't_coords', T, 'WriteMode', 'append');
        fprintf('\n');
    end
end
