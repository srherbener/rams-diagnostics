function [] = GenTotPcp(ConfigFile)
% GenTotPcp generate time series of total accumumlated precip

    % Read the config file to get the structure of how the data is laid out in
    % the file system.
    [ Config ] = ReadConfig(ConfigFile);

    Ddir = Config.DiagDir;
    Tdir = Config.TsavgDir;

    if (exist(Ddir, 'dir') ~= 7)
      mkdir(Ddir);
    end

    %PCPname = 'hda_totpcp';
    %PCPfprefix = 'hda_totpcp';
    PCPname = 'hda_accpr';
    PCPfprefix = 'hda_accpr';

    for icase = 1:length(Config.Cases)
        Case = Config.Cases(icase).Cname;
        PCPfile = sprintf('%s/%s_%s.h5', Tdir, PCPfprefix, Case);

        fprintf('***************************************************************\n');
        fprintf('Generating total precip data:\n');
        fprintf('  Case: %s\n', Case);
        fprintf('  Input precip file: %s\n', PCPfile);
        fprintf('    Var name: %s\n', PCPname);
        fprintf('\n');
 
        % precip will be organized as (x,t)
        %   x has size 2 and holds:
        %     x(1) --> sum of horizontal domain
        %     x(2) --> total number of points used to create y(1)
        %   t is time
        PCP = squeeze(hdf5read(PCPfile, PCPname));
        T   = squeeze(hdf5read(PCPfile, 't_coords'));

        % all we need is the time series in x(1) (sum across domain)
        TOT_PCP = squeeze(PCP(1,:));

        % Output coordinates
        X = 1;
        Y = 1;
        Z = 1;

        Nx = 1;
        Ny = 1;
        Nz = 1;
        Nt = length(T);

        % write out the series
        OutFile = sprintf('%s/totpcp_%s.h5', Ddir, Case);
        fprintf('Writing: %s\n', OutFile);

        Ovar = reshape(TOT_PCP, [ Nx Ny Nz Nt ]);

        hdf5write(OutFile, '/totpcp', Ovar);
        
        hdf5write(OutFile, 'x_coords', X, 'WriteMode', 'append');
        hdf5write(OutFile, 'y_coords', Y, 'WriteMode', 'append');
        hdf5write(OutFile, 'z_coords', Z, 'WriteMode', 'append');
        hdf5write(OutFile, 't_coords', T, 'WriteMode', 'append');
        fprintf('\n');
    end

end
