function [ ] = GenBuoyancyFlux(ConfigFile)
% GenBuoyancyFlux generate buoyancy from w' * theta_v' and theta-bar

    % Read the config file to get the structure of how the data is laid out in
    % the file system.
    [ Config ] = ReadConfig(ConfigFile);

    Ddir = Config.DiagDir;

    WpTvpName = 'w-theta_v';
    TbName = 'theta';

    WpTvpFprefix = 'w_thetav_M2_T12_T36';
    TbFprefix = 'theta_M1_T12_T36';

    g = 9.8;

    for icase = 1:length(Config.Cases)
        Case = Config.Cases(icase).Cname;

        WpTvpFile = sprintf('%s/%s_%s.h5', Ddir, WpTvpFprefix, Case);
        TbFile = sprintf('%s/%s_%s.h5', Ddir, TbFprefix, Case);
        OutFile = sprintf('%s/buoy_flux_T12_T36_%s.h5', Ddir, Case);

        fprintf('***************************************************************\n');
        fprintf('Generating buoyancy flux:\n');
        fprintf('  Case: %s\n', Case);
        fprintf('  Input virtual potential temp flux file: %s\n', WpTvpFile);
        fprintf('    Var name: %s\n', WpTvpName);
        fprintf('  Input mean potential temp file: %s\n', TbFile);
        fprintf('    Var name: %s\n', TbName);
        fprintf('\n');

        % Buoyancy flux:
        %   B = (g/ThetaBar) * (Wprime*ThetavPrime)
        %
        WP_TVP = squeeze(hdf5read(WpTvpFile, WpTvpName));
        TB     = squeeze(hdf5read(TbFile, TbName));

        % Grab coordinates: Z
        Z = hdf5read(WpTvpFile, 'z_coords');
        
        BF = (g ./ TB) .* WP_TVP;

        % output
        fprintf('Writing: %s\n', OutFile);
        hdf5write(OutFile, '/BUOY_FLUX', BF);
        hdf5write(OutFile, 'Z', Z, 'WriteMode', 'append');
        fprintf('\n');

    end

end
