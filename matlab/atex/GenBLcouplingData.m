function [ ] = GenBLcouplingData(ConfigFile)
% GenBLcouplingData generate time series of LCL and BL cloud base

    % Read the config file to get the structure of how the data is laid out in
    % the file system.
    [ Config ] = ReadConfig(ConfigFile);

    Ddir = Config.DiagDir;
    Tdir = Config.TsavgDir;
    Hdir = 'HDF5';

    UndefVal = Config.UndefVal;

    Tname = 'tempc';
    Tfprefix = 'sfc_tempc';
    TDname = 'dewptc';
    TDfprefix = 'sfc_dewptc';
    CldName = 'hda_cloud';
    CldFprefix = 'hda_cloud';

    for icase = 1:length(Config.Cases)
	Case = Config.Cases(icase).Cname;
	Tfile = sprintf('%s/%s-%s-AS-1999-02-10-040000-g1.h5', Hdir, Tfprefix, Case);
	TDfile = sprintf('%s/%s-%s-AS-1999-02-10-040000-g1.h5', Hdir, TDfprefix, Case);
	CldFile = sprintf('%s/%s_%s.h5', Tdir, CldFprefix, Case);

        fprintf('***************************************************************\n');
        fprintf('Generating BL coupling data:\n');
        fprintf('  Case: %s\n', Case);
        fprintf('  Input temp file: %s\n', Tfile);
        fprintf('    Var name: %s\n', Tname);
        fprintf('  Input dewpoint file: %s\n', TDfile);
        fprintf('    Var name: %s\n', TDname);
        fprintf('  Input cloud file: %s\n', CldFile);
        fprintf('    Var name: %s\n', CldName);
        fprintf('\n');
 
        % Cloud will be organized as (y,z,t)
        %   y has size 2 and holds:
        %     y(1) --> sum of horizontal domain
        %     y(2) --> total number of points used to create y(1)
        %   z is height
        %   t is time
        CLD = squeeze(hdf5read(CldFile, CldName));
        TEMP = squeeze(hdf5read(Tfile, Tname));
        TEMPD = squeeze(hdf5read(TDfile, TDname));

        % Grab height and time coordinates
        Z = squeeze(hdf5read(CldFile, 'z_coords'));
        T = squeeze(hdf5read(Tfile, 't_coords'));
        Nt = length(T);

        % Range to look for cloud base
        Z1 = find(Z >= 250, 1, 'first');
        Z2 = find(Z <= 2500, 1, 'last');
 
        % cloud base
        %   calculate average profile
        %   find the index of the first point (from sfc) that has cloud profile > 0.01 g/kg
        % CLD_PROF will be (z,t)
        CLD_PROF = squeeze(CLD(1,Z1:Z2,:) ./ CLD(2,Z1:Z2,:));
        CLD_BASE = zeros([ 1 Nt ]);
        for it = 1:Nt
          ZB = find(squeeze(CLD_PROF(:,it) > 0.01), 1, 'first');
          if (isempty(ZB))
            CLD_BASE(it) = nan;
          else
            ZB = ZB + Z1 - 1;
            CLD_BASE(it) = Z(ZB);
          end
        end

        % average LCL
        T_DIFF = TEMP - TEMPD;
        LCL = 125 .* squeeze(mean(mean(T_DIFF,1),2));

        % output --> Use REVU format, 4D var, *_coords
        X = 1;
        Y = 1;
        Z = 1;

        OutFile = sprintf('%s/BL_coupling_%s.h5', Ddir, Case);
        fprintf('Writing: %s\n', OutFile);
 
        Ovar = reshape(LCL, [ 1 1 1 Nt ]);
        hdf5write(OutFile, '/avg_lcl', Ovar);

        Ovar = reshape(CLD_BASE, [ 1 1 1 Nt ]);
        hdf5write(OutFile, '/avg_cloud_base', Ovar, 'WriteMode', 'append');

        hdf5write(OutFile, 'x_coords', X, 'WriteMode', 'append');
        hdf5write(OutFile, 'y_coords', Y, 'WriteMode', 'append');
        hdf5write(OutFile, 'z_coords', Z, 'WriteMode', 'append');
        hdf5write(OutFile, 't_coords', T, 'WriteMode', 'append');
    fprintf('\n');
    end
end
