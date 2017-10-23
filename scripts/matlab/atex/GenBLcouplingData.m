function [ ] = GenBLcouplingData(ConfigFile)
% GenBLcouplingData generate time series of LCL and BL cloud base

    % Mixing ratio for determining cloud base
    CLD_THRESHOLD = 0.1;   % g/kg
    Nfilter = 5;           % use five point in running average smoothing filter

    CBASE_THRESHOLD = 0.5; % 50% of max value in PDF

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
    CldName = 'hist_cloud';
    CldFprefix = 'hist_cloud_m2';

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
 
        % Use nctoolbox in order to walk through the data per time step
        CLD_DS = ncgeodataset(CldFile);
        CLD_VAR = CLD_DS.geovariable(CldName);
        X_VAR = CLD_DS.geovariable('x_coords');
        Z_VAR = CLD_DS.geovariable('z_coords');
        T_VAR = CLD_DS.geovariable('t_coords');

        X = X_VAR.data(:);
        Z = Z_VAR.data(:);
        T = T_VAR.data(:);

        Nt = length(T);

        TEMP_DS = ncgeodataset(Tfile);
        TEMP_VAR = TEMP_DS.geovariable(Tname);

        TEMPD_DS = ncgeodataset(TDfile);
        TEMPD_VAR = TEMPD_DS.geovariable(TDname);

        % Find cloud base of the stratiform deck
        %   Use histogram of cloud mixing ratio
        %   Sum up the counts from the bins >= CLD_THRESHOLD which
        %     will create a vertical profile of counts of cloudy regions
        %   Smooth the profile
        %   Change profile into PDF
        %   Select first occurance of PDF entry, starting from ground, that crosses CBASE_THRESHOLD
        %
        % Note that out of the files, the variables will have the form (t, z, y, x)
        %     CLD will have dummy variable (size == 1) y  --> (t, z, 1, x)
        %     TEMP and TEMPD will have dummy variable z   --> (t, 1, y, x)
        CLD_BASE = zeros([ 1 Nt ]);
        LCL = zeros([ 1 Nt ]);

        X1 = find(X >= CLD_THRESHOLD, 1, 'first');

        for it = 1:Nt
          % Can get the case where all counts are zero, in this case CLD_PDF will
          % come out all nans. For this case place a nan in CLD_BASE.
          CLD = squeeze(CLD_VAR.data(it, :, :, :));       % CLD will be (z,x)
          CLD_PDF = smooth(squeeze(sum(CLD(:,X1:end), 2)), Nfilter);
          SUM = sum(CLD_PDF);
          if (SUM == 0)
            CLD_BASE(it) = nan;
          else
            CLD_PDF = CLD_PDF ./ SUM;
            ZB = find(CLD_PDF >= max(CLD_PDF)*CBASE_THRESHOLD, 1, 'first');
            CLD_BASE(it) = Z(ZB);
          end
        
          % average LCL
          TEMP  = TEMP_VAR.data(it,:,:,:);    % TEMP, TEMPD will be (y,x)
          TEMPD = TEMPD_VAR.data(it,:,:,:);
          Tmean = mean(TEMP(:))   + 273.15; % convert to Kelvin
          TDmean = mean(TEMPD(:)) + 273.15;

          % temp of adiabatically lifted parcel at condensation level
          % formula is from Bolton, 1980 MWR "The Computation of Equivalent Potential Temperature"
          %    equation (15)
          Tlcl = 1 / ((1/(TDmean-56)) + (log(Tmean/TDmean)/800)) + 56;

          % the change in height will just be (T - Tlcl) divided
          % by the adiabatic lapse rate: 9.8 K / km -> 0.0098 K / m -> 102 m / K
          LCL(it) = Z(2) + ((Tmean - Tlcl) * 102);
        end

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
