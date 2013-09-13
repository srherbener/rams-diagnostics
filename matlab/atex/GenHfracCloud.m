function [ ] = GenHfracCloud(ConfigFile)
% GenHfracCloud generate horizontal domain cloud fraction (profiles)

    % Read the config file to get the structure of how the data is laid out in
    % the file system.
    [ Config ] = ReadConfig(ConfigFile);

    Ddir = Config.DiagDir;
    Tdir = Config.TsavgDir;

    CfracName = 'hfrac_cloud';
    CfracFprefix = 'hfrac_cloud';

    Time1 = 12;
    Time2 = 36;
    TimeStr = sprintf('T%d_T%d', Time1, Time2);

    for icase = 1:length(Config.Cases)
        Case = Config.Cases(icase).Cname;

        CfracFile = sprintf('%s/%s_%s.h5', Tdir, CfracFprefix, Case);
        OutFile = sprintf('%s/hfrac_cloud_%s_%s.h5', Ddir, TimeStr, Case);

        fprintf('***************************************************************\n');
        fprintf('Generating cloud fraction profiles:\n');
        fprintf('  Case: %s\n', Case);
        fprintf('  Input cloud fraction file: %s\n', CfracFile);
        fprintf('    Var name: %s\n', CfracName);
        fprintf('\n');

        % Data will be (z,t) after squeeze
        CF = squeeze(hdf5read(CfracFile, CfracName));

        % Grab height coordinates
        Z = hdf5read(CfracFile, 'z_coords');
        Nz = length(Z);
        
        % Take time average
        Time = hdf5read(CfracFile, 't_coords')/3600; % hrs
        T1 = find(Time >= Time1, 1, 'first');
        T2 = find(Time <= Time2, 1, 'last');
        CF_AVG = squeeze(mean(CF(:,T1:T2),2));

        % output --> Use REVU format, 4D var, *_coords
        % fabricate x, y, t coords
        X = 1;
        Y = 1;
        T = 1;
        Ovar = reshape(CF_AVG, [ 1 1 Nz 1 ]);
        
        fprintf('Writing: %s\n', OutFile);
        hdf5write(OutFile, '/hfrac_cloud', Ovar);
        hdf5write(OutFile, 'x_coords', X, 'WriteMode', 'append');
        hdf5write(OutFile, 'y_coords', Y, 'WriteMode', 'append');
        hdf5write(OutFile, 'z_coords', Z, 'WriteMode', 'append');
        hdf5write(OutFile, 't_coords', T, 'WriteMode', 'append');
        fprintf('\n');
    end
end
