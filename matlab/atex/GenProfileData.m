function [ ] = GenProfileData(ConfigFile)
% GenProfileData generate data indexed by SST and CCN concentration for
% profiles

    % Read the config file to get the structure of how the data is laid out in
    % the file system.
    [ Config ] = ReadConfig(ConfigFile);

    Ddir = Config.DiagDir;
    Tdir = Config.TsavgDir;

    UndefVal = Config.UndefVal;

    CldName = 'hda_cloud';
    CldFprefix = 'hda_cloud';
    CldFprefix_lwp_0p01 = 'hda_cloud_lwp_0p01';
    CldFprefix_lwp_0p10 = 'hda_cloud_lwp_0p10';
    CldFprefix_lwp_1p00 = 'hda_cloud_lwp_1p00';

    Tstart = 12;
    Tend = 36;

    i = 0;
    for icase = 1:length(Config.Cases)
        Case = Config.Cases(icase).Cname;
        CldFile = sprintf('%s/%s_%s.h5', Tdir, CldFprefix, Case);
        CldFile_lwp_0p01 = sprintf('%s/%s_%s.h5', Tdir, CldFprefix_lwp_0p01, Case);
        CldFile_lwp_0p10 = sprintf('%s/%s_%s.h5', Tdir, CldFprefix_lwp_0p10, Case);
        CldFile_lwp_1p00 = sprintf('%s/%s_%s.h5', Tdir, CldFprefix_lwp_1p00, Case);

        fprintf('***************************************************************\n');
        fprintf('Generating profile data:\n');
        fprintf('  Case: %s\n', Case);
        fprintf('  Input cloud file: %s\n', CldFile);
        fprintf('  Input cloud file: %s\n', CldFile_lwp_0p01);
        fprintf('  Input cloud file: %s\n', CldFile_lwp_0p10);
        fprintf('  Input cloud file: %s\n', CldFile_lwp_1p00);
        fprintf('    Var name: %s\n', CldName);
        fprintf('  Start Time: %d\n', Tstart);
        fprintf('  End Time: %d\n', Tend);
        fprintf('\n');
 
        % Cloud will be organized as (z,t)
        CLD = squeeze(hdf5read(CldFile, CldName));
        CLD_LWP_0P01 = squeeze(hdf5read(CldFile_lwp_0p01, CldName));
        CLD_LWP_0P10 = squeeze(hdf5read(CldFile_lwp_0p10, CldName));
        CLD_LWP_1P00 = squeeze(hdf5read(CldFile_lwp_1p00, CldName));

        % Grab time coordinates
        T = hdf5read(CldFile, 't_coords')/3600; % hrs
 
        % Select out the given time step
        T1 = find(T >= Tstart, 1, 'first');
        T2 = find(T <= Tend, 1, 'last');

        % Record data
        i = i + 1;

        % compute temporal average profile
        % put in nan's for the undefined value and use nanmean
        CLD(CLD == UndefVal) = nan;
        CLD_LWP_0P01(CLD_LWP_0P01 == UndefVal) = nan;
        CLD_LWP_0P10(CLD_LWP_0P10 == UndefVal) = nan;
        CLD_LWP_1P00(CLD_LWP_1P00 == UndefVal) = nan;

        CLD_PROF = nanmean(CLD(:,T1:T2),2);
        CLD_PROF_LWP_0P01 = nanmean(CLD_LWP_0P01(:,T1:T2),2);
        CLD_PROF_LWP_0P10 = nanmean(CLD_LWP_0P10(:,T1:T2),2);
        CLD_PROF_LWP_1P00 = nanmean(CLD_LWP_1P00(:,T1:T2),2);
        
            % output --> Use REVU format, 4D var, *_coords
        X = 1;
        Y = 1;
        Z = squeeze(hdf5read(CldFile, 'z_coords'));
        T = 1;
        
        Nz = length(Z);
        Ovar = reshape(CLD_PROF, [ 1 1 Nz 1 ]);
        Ovar_lwp_0p01 = reshape(CLD_PROF_LWP_0P01, [ 1 1 Nz 1 ]);
        Ovar_lwp_0p10 = reshape(CLD_PROF_LWP_0P10, [ 1 1 Nz 1 ]);
        Ovar_lwp_1p00 = reshape(CLD_PROF_LWP_1P00, [ 1 1 Nz 1 ]);

        OutFile = sprintf('%s/cloud_profile_T%d_T%d_%s.h5', Ddir, Tstart, Tend, Case);
        fprintf('Writing: %s\n', OutFile);
        hdf5write(OutFile, '/cloud', Ovar);
        hdf5write(OutFile, '/cloud_lwp_0p01', Ovar_lwp_0p01, 'WriteMode', 'append');
        hdf5write(OutFile, '/cloud_lwp_0p10', Ovar_lwp_0p10, 'WriteMode', 'append');
        hdf5write(OutFile, '/cloud_lwp_1p00', Ovar_lwp_1p00, 'WriteMode', 'append');
        
        hdf5write(OutFile, 'x_coords', X, 'WriteMode', 'append');
        hdf5write(OutFile, 'y_coords', Y, 'WriteMode', 'append');
        hdf5write(OutFile, 'z_coords', Z, 'WriteMode', 'append');
        hdf5write(OutFile, 't_coords', T, 'WriteMode', 'append');
    fprintf('\n');
    end
end
