function [ ] = GenAlbedo(ConfigFile)
% GenAlbedo generate cloud albedo from cloud thickness

    % Read the config file to get the structure of how the data is laid out in
    % the file system.
    [ Config ] = ReadConfig(ConfigFile);

    Hdir = 'HDF5';
    Ddir = Config.DiagDir;

    CtName = 'cloud_opt_thick';
    CtFprefix = 'cloud_opt_thick';

    Times = [ 12 18 24 30 36 ];

    g = 0.85; % scattering asymmetry parameter (cloud droplets interacting with visible light)

    for icase = 1:length(Config.Cases)
        Case = Config.Cases(icase).Cname;
        CtFile = sprintf('%s/%s-%s-AS-1999-02-10-040000-g1.h5', Hdir, CtFprefix, Case);

        for it = 1:length(Times)
            Time = Times(it);
            TimeSec = Time * 3600;

            OutFile = sprintf('%s/albedo_T%d_%s.h5', Ddir, Time, Case);
    
            fprintf('***************************************************************\n');
            fprintf('Generating albedo:\n');
            fprintf('  Case: %s\n', Case);
            fprintf('  Input cloud optical thickness file: %s\n', CtFile);
            fprintf('    Var name: %s\n', CtName);
            fprintf('  Time: %d\n', Time);
            fprintf('\n');
    
           % Albedo:
           %   A = (1-g)*CT / (((1-g)*CT) + 2)
           %
           HDATA = squeeze(hdf5read(CtFile, CtName));
   
           % Grab coordinates
           X = hdf5read(CtFile, 'x_coords');
           Y = hdf5read(CtFile, 'y_coords');
           Z = hdf5read(CtFile, 'z_coords');
           T = hdf5read(CtFile, 't_coords');
           Nx = length(X);
           Ny = length(Y);
           Nz = length(Z);
           Nt = 1;         % after filtering

           % Select out the given time step
           % HDATA will be (x,y,t)
           T1 = find(T >= TimeSec, 1, 'first');
           CT = squeeze(HDATA(:,:,T1));
           
           % albedo formula from Borhen, 1987
           ALB = ((1 - g) .* CT) ./ (((1 - g) .* CT) + 2);
   
           % output --> Use REVU format, 4D var, *_coords
           % fabricate x, y, t coords
           Ovar = reshape(ALB, [ Nx Ny Nz Nt ]);
            
           fprintf('Writing: %s\n', OutFile);
           hdf5write(OutFile, '/albedo', Ovar);
           hdf5write(OutFile, 'x_coords', X, 'WriteMode', 'append');
           hdf5write(OutFile, 'y_coords', Y, 'WriteMode', 'append');
           hdf5write(OutFile, 'z_coords', Z, 'WriteMode', 'append');
           hdf5write(OutFile, 't_coords', TimeSec, 'WriteMode', 'append');
           fprintf('\n');
        end
    end
end
