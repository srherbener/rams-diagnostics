function [ ] = GenSfcDiv(ConfigFile)
% GenSfcDif generate surface divergence

    % Read the config file to get the structure of how the data is laid out in
    % the file system.
    [ Config ] = ReadConfig(ConfigFile);

    Hdir = 'HDF5';
    Ddir = Config.DiagDir;

    SdName = 'sfc_div';
    SdFprefix = 'sfc_div';

    Times = [ 12 18 24 30 36 ];

    g = 0.85; % scattering asymmetry parameter (cloud droplets interacting with visible light)

    for icase = 1:length(Config.Cases)
        Case = Config.Cases(icase).Cname;
        SdFile = sprintf('%s/%s-%s-AS-1999-02-10-040000-g1.h5', Hdir, SdFprefix, Case);

        for it = 1:length(Times)
            Time = Times(it);
            TimeSec = Time * 3600;

            OutFile = sprintf('%s/sfc_div_T%d_%s.h5', Ddir, Time, Case);
    
            fprintf('***************************************************************\n');
            fprintf('Generating surface divergence:\n');
            fprintf('  Case: %s\n', Case);
            fprintf('  Input surface divergence file: %s\n', SdFile);
            fprintf('    Var name: %s\n', SdName);
            fprintf('  Time: %d\n', Time);
            fprintf('\n');
    
           % simply select the data and dump into the output file
           HDATA = squeeze(hdf5read(SdFile, SdName));
   
           % Grab coordinates
           X = hdf5read(SdFile, 'x_coords');
           Y = hdf5read(SdFile, 'y_coords');
           Z = hdf5read(SdFile, 'z_coords');
           T = hdf5read(SdFile, 't_coords');
           Nx = length(X);
           Ny = length(Y);
           Nz = length(Z);
           Nt = 1;         % after filtering

           % Select out the given time step
           % HDATA will be (x,y,t)
           T1 = find(T >= TimeSec, 1, 'first');
           SDIV = squeeze(HDATA(:,:,T1));
           
           % output --> Use REVU format, 4D var, *_coords
           % fabricate x, y, t coords
           Ovar = reshape(SDIV, [ Nx Ny Nz Nt ]);
            
           fprintf('Writing: %s\n', OutFile);
           hdf5write(OutFile, '/sfc_div', Ovar);
           hdf5write(OutFile, 'x_coords', X, 'WriteMode', 'append');
           hdf5write(OutFile, 'y_coords', Y, 'WriteMode', 'append');
           hdf5write(OutFile, 'z_coords', Z, 'WriteMode', 'append');
           hdf5write(OutFile, 't_coords', TimeSec, 'WriteMode', 'append');
           fprintf('\n');
        end
    end
end
