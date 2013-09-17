function [ ] = GenProfileData(ConfigFile)
% GenProfileData generate data indexed by SST and CCN concentration for
% profiles

    % Read the config file to get the structure of how the data is laid out in
    % the file system.
    [ Config ] = ReadConfig(ConfigFile);

    Ddir = Config.DiagDir;
    Tdir = Config.TsavgDir;

    CldName = 'hda_cloud';
    CldFprefix = 'hda_cloud';

    Tstart = 12;
    Tend = 36;

    i = 0;
    for icase = 1:length(Config.Cases)
        Case = Config.Cases(icase).Cname;
        CldFile = sprintf('%s/%s_%s.h5', Tdir, CldFprefix, Case);

        fprintf('***************************************************************\n');
        fprintf('Generating profile data:\n');
        fprintf('  Case: %s\n', Case);
        fprintf('  Input cloud file: %s\n', CldFile);
        fprintf('    Var name: %s\n', CldName);
        fprintf('  Start Time: %d\n', Tstart);
        fprintf('  End Time: %d\n', Tend);
        fprintf('\n');
 
        % Cloud will be organized as (z,t)
        CLD = squeeze(hdf5read(CldFile, CldName));
    
        % Grab time coordinates
        T = hdf5read(CldFile, 't_coords')/3600; % hrs
 
        % Select out the given time step
        T1 = find(T >= Tstart, 1, 'first');
        T2 = find(T <= Tend, 1, 'last');

        % Record data
        i = i + 1;

        % compute temporal average profile
        CLD_PROF = mean(CLD(:,T1:T2),2);
        
            % output --> Use REVU format, 4D var, *_coords
        X = 1;
        Y = 1;
        Z = squeeze(hdf5read(CldFile, 'z_coords'));
        T = 1;
        
        Nz = length(Z);
        Ovar = reshape(CLD_PROF, [ 1 1 Nz 1 ]);

        OutFile = sprintf('%s/cloud_profile_T%d_T%d_%s.h5', Ddir, Tstart, Tend, Case);
        fprintf('Writing: %s\n', OutFile);
        hdf5write(OutFile, '/cloud', Ovar);
        
        hdf5write(OutFile, 'x_coords', X, 'WriteMode', 'append');
        hdf5write(OutFile, 'y_coords', Y, 'WriteMode', 'append');
        hdf5write(OutFile, 'z_coords', Z, 'WriteMode', 'append');
        hdf5write(OutFile, 't_coords', T, 'WriteMode', 'append');
    fprintf('\n');
    end
end
