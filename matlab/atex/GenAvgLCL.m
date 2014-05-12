function [ ] = GenAvgLCL(ConfigFile)
% GenAvgLCL generate time series of domain average LCL

    % Read the config file to get the structure of how the data is laid out in
    % the file system.
    [ Config ] = ReadConfig(ConfigFile);

    % level where sfc_temp and sfc_dewptc are measured
    Z_BASE = 50;  % 50m AGL

    BINS_LCL = 0:10:4000;
    Nb = length(BINS_LCL);

    Ddir = Config.DiagDir;
    Tdir = Config.TsavgDir;
    Fdir = 'FILTERS';
    Hdir = 'HDF5';

    UndefVal = Config.UndefVal;

    Tname = 'tempc';
    Tfprefix = 'sfc_tempc';
    TDname = 'dewptc';
    TDfprefix = 'sfc_dewptc';

    Fname    = 'filter';
    Ffprefix = 'stall';

    for icase = 1:length(Config.Cases)
	Case = Config.Cases(icase).Cname;
	Tfile = sprintf('%s/%s-%s-AS-1999-02-10-040000-g1.h5', Hdir, Tfprefix, Case);
	TDfile = sprintf('%s/%s-%s-AS-1999-02-10-040000-g1.h5', Hdir, TDfprefix, Case);
	Ffile = sprintf('%s/%s_%s.h5', Fdir, Ffprefix, Case);

        fprintf('***************************************************************\n');
        fprintf('Generating horizontal domain average LCL:\n');
        fprintf('  Case: %s\n', Case);
        fprintf('  Input temp file: %s\n', Tfile);
        fprintf('    Var name: %s\n', Tname);
        fprintf('  Input dewpoint file: %s\n', TDfile);
        fprintf('    Var name: %s\n', TDname);
        fprintf('  Input filter file: %s\n', Ffile);
        fprintf('    Var name: %s\n', Fname);
        fprintf('\n');


        TEMP_ALL   = squeeze(hdf5read(Tfile, Tname));
        TEMPD_ALL  = squeeze(hdf5read(TDfile, TDname));
        FILTER_ALL = squeeze(hdf5read(Ffile, Fname));

        T = squeeze(hdf5read(Tfile, 't_coords'));
        Nt = length(T);

        % Mimic the manner in which the "hda" averaging function in tsavg creates
        % its data so that the output of this function can be processed in the same
        % manner as tsavg hda output.
        %
        % Output needs to be 4D (x,y,z,t) where x and z have size 1, y has size 2
        % and t has size Nt. y(1) holds the sums and y(2) holds the number of points.
        HDA_LCL  = zeros([ 2 Nt ]);
        HIST_LCL  = zeros([ Nb Nt ]);
        for it = 1:Nt
          % record the sum of LCL values across the domain
          % and the number of points so that averages can
          % formed at each time step or across multiple
          % time steps
          %
          % exclude the lateral boundaries
          TEMP   = squeeze(TEMP_ALL(2:end-1,2:end-1,it)) + 273.15;   % TEMP, TEMPD will be (x,y)
          TEMPD  = squeeze(TEMPD_ALL(2:end-1,2:end-1,it)) + 273.15;  % convert to Kelvin
          FILTER = squeeze(FILTER_ALL(2:end-1,2:end-1,it));

          % temp of adiabatically lifted parcel at condensation level
          % formula is from Bolton, 1980 MWR "The Computation of Equivalent Potential Temperature"
          %    equation (15)
          T_LCL = 1 ./ ((1./(TEMPD-56)) + (log(TEMP./TEMPD)./800)) + 56;

          % the change in height will just be (T - Tlcl) divided
          % by the adiabatic lapse rate: 9.8 K / km -> 0.0098 K / m -> 102 m / K
          LCL = Z_BASE + ((TEMP - T_LCL) .* 102);

          % record sum and npts
          LCL = LCL(:);
          FILTER = FILTER(:);

          LCL_FILT = LCL(FILTER == 1);
   
          HDA_LCL(1,it) = sum(LCL);
          HDA_LCL(2,it) = length(LCL);
          HIST_LCL(:,it) = histc(LCL, BINS_LCL);

          HDA_LCL_FILT(1,it) = sum(LCL_FILT);
          HDA_LCL_FILT(2,it) = length(LCL_FILT);
          HIST_LCL_FILT(:,it) = histc(LCL_FILT, BINS_LCL);
        end

        % output --> Use REVU format, 4D var, *_coords
        Y = [ 1 2 ];
        Z = 50; % first model level above ground

        OutFile = sprintf('%s/hda_lcl_%s.h5', Ddir, Case);
        fprintf('Writing: %s\n', OutFile);
 
        Ovar = reshape(HDA_LCL, [ 1 2 1 Nt ]);
        hdf5write(OutFile, '/hda_lcl', Ovar);

        Ovar = reshape(HIST_LCL, [ Nb 1 1 Nt ]);
        hdf5write(OutFile, '/hist_lcl', Ovar, 'WriteMode', 'append');

        hdf5write(OutFile, 'x_coords', BINS_LCL, 'WriteMode', 'append');
        hdf5write(OutFile, 'y_coords', Y, 'WriteMode', 'append');
        hdf5write(OutFile, 'z_coords', Z, 'WriteMode', 'append');
        hdf5write(OutFile, 't_coords', T, 'WriteMode', 'append');
        fprintf('\n');


        OutFile = sprintf('%s/hda_lcl_stall_%s.h5', Ddir, Case);
        fprintf('Writing: %s\n', OutFile);

        Ovar = reshape(HDA_LCL_FILT, [ 1 2 1 Nt ]);
        hdf5write(OutFile, '/hda_lcl', Ovar);

        Ovar = reshape(HIST_LCL_FILT, [ Nb 1 1 Nt ]);
        hdf5write(OutFile, '/hist_lcl', Ovar, 'WriteMode', 'append');

        hdf5write(OutFile, 'x_coords', BINS_LCL, 'WriteMode', 'append');
        hdf5write(OutFile, 'y_coords', Y, 'WriteMode', 'append');
        hdf5write(OutFile, 'z_coords', Z, 'WriteMode', 'append');
        hdf5write(OutFile, 't_coords', T, 'WriteMode', 'append');
        fprintf('\n');
    end
end
