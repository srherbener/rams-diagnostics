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

    Filters = {
      'none'
      'lwp0p01'
      'lwp0p10'
      'lwp1p00'
      };

    g = 9.8;

    for icase = 1:length(Config.Cases)
        Case = Config.Cases(icase).Cname;

        for ifilt = 1:length(Filters)
            Filter = Filters{ifilt};

            if (strcmp(Filter, 'none'))
                WpTvpFile = sprintf('%s/%s_%s.h5', Ddir, WpTvpFprefix, Case);
                TbFile = sprintf('%s/%s_%s.h5', Ddir, TbFprefix, Case);
                OutFile = sprintf('%s/buoy_flux_T12_T36_%s.h5', Ddir, Case);
            else
                WpTvpFile = sprintf('%s/%s_%s_%s.h5', Ddir, WpTvpFprefix, Filter, Case);
                TbFile = sprintf('%s/%s_%s_%s.h5', Ddir, TbFprefix, Filter, Case);
                OutFile = sprintf('%s/buoy_flux_%s_T12_T36_%s.h5', Ddir, Filter, Case);
            end
    
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
    
            % Grab height coordinates
            Z = hdf5read(WpTvpFile, 'z_coords');
            Nz = length(Z);
            
            BF = (g ./ TB) .* WP_TVP;
    
            % output --> Use REVU format, 4D var, *_coords
            % fabricate x, y, t coords
            X = 1;
            Y = 1;
            T = 1;
            Ovar = reshape(BF, [ 1 1 Nz 1 ]);
            
            fprintf('Writing: %s\n', OutFile);
            hdf5write(OutFile, '/BUOY_FLUX', Ovar);
            hdf5write(OutFile, 'x_coords', X, 'WriteMode', 'append');
            hdf5write(OutFile, 'y_coords', Y, 'WriteMode', 'append');
            hdf5write(OutFile, 'z_coords', Z, 'WriteMode', 'append');
            hdf5write(OutFile, 't_coords', T, 'WriteMode', 'append');
            fprintf('\n');
        end
    end
end
