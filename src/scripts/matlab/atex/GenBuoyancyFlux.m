function [ ] = GenBuoyancyFlux(ConfigFile)
% GenBuoyancyFlux generate buoyancy from w' * theta_v' and theta-bar

    % Read the config file to get the structure of how the data is laid out in
    % the file system.
    [ Config ] = ReadConfig(ConfigFile);

    Ddir = Config.DiagDir;

    VarSets = {
      { 'moments' { 'col_up0p10_all' 'col_dn0p10_all' 'col_ud0p10_all' 'all_cld_all' 'up0p10_all_cld_all' 'dn0p10_all_cld_all' 'ud0p10_all_cld_all' 'up0p01_all_cld_all' 'dn0p01_all_cld_all' 'ud0p01_all_cld_all' } 'buoy_flux' }
      };

    Nvsets = length(VarSets);

    g = 9.8;

    for icase = 1:length(Config.Cases)
        Case = Config.Cases(icase).Cname;

        for ivset = 1:Nvsets
            InFilePrefix  = VarSets{ivset}{1};
            InVlist       = VarSets{ivset}{2};
            OutFilePrefix = VarSets{ivset}{3};

            InFile = sprintf('%s/%s_%s.h5', Ddir, InFilePrefix, Case);
            OutFile = sprintf('%s/%s_%s.h5', Ddir, OutFilePrefix, Case);

            fprintf('***************************************************************\n');
            fprintf('Generating buoyancy flux:\n');
            fprintf('  Case: %s\n', Case);
            fprintf('  Reading: %s\n', InFile);
            fprintf('  Writing: %s\n', OutFile);
            fprintf('\n');

            % Write a header so that a new file is created, and 
            % the writes in the following loop can use 'append' mode.
            hdf5write(OutFile, '/Header', 'Buoyancy Flux');

            for ivar = 1:length(InVlist)
              InVsuffix = InVlist{ivar};

              TbVname    = sprintf('theta_v_%s', InVsuffix);   % mean theta_v
              TvpVname   = sprintf('w-theta_v_%s', InVsuffix); % w-theta_v covariance
              BfluxVname = sprintf('buoy_flux_%s', InVsuffix); % buoyance flux
    
              fprintf('    Vars: %s, %s --> %s\n', TvpVname, TbVname, BfluxVname);
    
              % Buoyancy flux:
              %   B = (g/ThetaBar) * (Wprime*ThetavPrime)
              %
              WP_TVP = squeeze(hdf5read(InFile, TvpVname));
              TB     = squeeze(hdf5read(InFile, TbVname));
      
              % Grab height coordinates
              Z = hdf5read(InFile, 'z_coords');
              Nz = length(Z);
              
              BF = (g ./ TB) .* WP_TVP;
      
              % output --> Use REVU format, 4D var, *_coords
              Ovar = reshape(BF, [ 1 1 Nz 1 ]);
              
              hdf5write(OutFile, BfluxVname, Ovar, 'WriteMode', 'append');
            end

            % write out the coordinates
            % fabricate x, y, t coords
            X = 1;
            Y = 1;
            T = 1;

            hdf5write(OutFile, 'x_coords', X,    'WriteMode', 'append');
            hdf5write(OutFile, 'y_coords', Y,    'WriteMode', 'append');
            hdf5write(OutFile, 'z_coords', Z,    'WriteMode', 'append');
            hdf5write(OutFile, 't_coords', T,    'WriteMode', 'append');
            fprintf('\n');
        end
    end
end
