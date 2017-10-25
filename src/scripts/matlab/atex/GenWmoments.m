function [ ] = GenWmoments(ConfigFile)
% GenWmoments generate w moments

    % Read the config file to get the structure of how the data is laid out in
    % the file system.
    [ Config ] = ReadConfig(ConfigFile);

    Ddir = Config.DiagDir;

    W1name = 'w';
    W2name = 'w-w';
    W3name = 'w-w-w';

    W1fprefix = 'w_M1_T12_T36';
    W2fprefix = 'w_M2_T12_T36';
    W3fprefix = 'w_M3_T12_T36';

    Filters = {
      'none'
      'lwp0p01'
      'lwp0p10'
      'lwp1p00'
      };

    for icase = 1:length(Config.Cases)
        Case = Config.Cases(icase).Cname;

        for ifilt = 1:length(Filters)
            Filter = Filters{ifilt};

            if (strcmp(Filter, 'none'))
                W1file = sprintf('%s/%s_%s.h5', Ddir, W1fprefix, Case);
                W2file = sprintf('%s/%s_%s.h5', Ddir, W2fprefix, Case);
                W3file = sprintf('%s/%s_%s.h5', Ddir, W3fprefix, Case);
                OutFile = sprintf('%s/w_moments_T12_T36_%s.h5', Ddir, Case);
            else
                W1file = sprintf('%s/%s_%s_%s.h5', Ddir, W1fprefix, Filter, Case);
                W2file = sprintf('%s/%s_%s_%s.h5', Ddir, W2fprefix, Filter, Case);
                W3file = sprintf('%s/%s_%s_%s.h5', Ddir, W3fprefix, Filter, Case);
                OutFile = sprintf('%s/w_moments_%s_T12_T36_%s.h5', Ddir, Filter, Case);
            end
    
            fprintf('***************************************************************\n');
            fprintf('Generating buoyancy flux:\n');
            fprintf('  Case: %s\n', Case);
            fprintf('  Input w first moment file: %s\n', W1file);
            fprintf('    Var name: %s\n', W1name);
            fprintf('  Input w second moment file: %s\n', W2file);
            fprintf('    Var name: %s\n', W2name);
            fprintf('  Input w third moment file: %s\n', W3file);
            fprintf('    Var name: %s\n', W3name);
            fprintf('\n');
    
            W1 = squeeze(hdf5read(W1file, W1name));
            W2 = squeeze(hdf5read(W2file, W2name));
            W3 = squeeze(hdf5read(W3file, W3name));
    
            % Grab height coordinates
            Z = hdf5read(W1file, 'z_coords');
            Nz = length(Z);
            
            % output --> Use REVU format, 4D var, *_coords
            % fabricate x, y, t coords
            X = 1;
            Y = 1;
            T = 1;
            W1out = reshape(W1, [ 1 1 Nz 1 ]);
            W2out = reshape(W2, [ 1 1 Nz 1 ]);
            W3out = reshape(W3, [ 1 1 Nz 1 ]);
            
            fprintf('Writing: %s\n', OutFile);
            hdf5write(OutFile, '/W1', W1out);
            hdf5write(OutFile, '/W2', W2out, 'WriteMode', 'append');
            hdf5write(OutFile, '/W3', W3out, 'WriteMode', 'append');
            hdf5write(OutFile, 'x_coords', X, 'WriteMode', 'append');
            hdf5write(OutFile, 'y_coords', Y, 'WriteMode', 'append');
            hdf5write(OutFile, 'z_coords', Z, 'WriteMode', 'append');
            hdf5write(OutFile, 't_coords', T, 'WriteMode', 'append');
            fprintf('\n');
        end
    end
end
