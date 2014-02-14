function [ ] = GenMomentData(ConfigFile)
% GenMomentData generate w moment and w flux data

    % Read the config file to get the structure of how the data is laid out in
    % the file system.
    [ Config ] = ReadConfig(ConfigFile);

    Ddir = Config.DiagDir;

    FprefixList = {
      'cloud_M1'
      'cloud_M1_c0p01'
      'cloud_M1_c0p10'
      'theta_e_M1'
      'w_M3'
      'w_speed_flux'
      'w_theta_flux'
      'w_theta_v_flux'
      'w_vapor_flux'
      };

    VnameList = {
      'cloud'
      'cloud'
      'cloud'
      'theta_e'
      'w-w-w'
      'w-speed'
      'w-theta'
      'w-theta_v'
      'w-vapor'
      };

    Tstart = 12;
    Tend = 36;

    for icase = 1:length(Config.Cases)
      Case = Config.Cases(icase).Cname;

      for ivar = 1:length(FprefixList)
        Fprefix = FprefixList{ivar};
        Fname = sprintf('%s/%s_%s.h5', Ddir, Fprefix, Case);
        Vname = VnameList{ivar};

        fprintf('***************************************************************\n');
        fprintf('Generating moment/flux profiles:\n');
        fprintf('  Case: %s\n', Case);
        fprintf('  Input file: %s\n', Fname);
        fprintf('    Var name: %s\n', Vname);
        fprintf('\n');

        TERMS = squeeze(hdf5read(Fname, Vname));
        NPTS = squeeze(hdf5read(Fname, 'num_points'));
        Z = squeeze(hdf5read(Fname, 'z_coords'));
        T = squeeze(hdf5read(Fname, 't_coords')) / 3600;   % hr
        
        Nz = length(Z);
        Nt = length(T);

        % *_M1 files will lose their first two dimensions, so put them back
        if (strcmp(Vname, 'cloud') || strcmp(Vname, 'theta_e'))
          TERMS = reshape(TERMS, [ 1 1 Nz Nt ]);
        end

        % select and sum up across the time range
        T1 = find(T >= Tstart, 1, 'first');
        T2 = find(T <= Tend, 1, 'last');

        [ MOMENTS OUT_NPTS ] = GenMoments(TERMS, NPTS, T1, T2);

        % MOMENTS is organized: (Nterms, Norder, Nz)
        %
        % For cloud and theta_e, MOMENTS will be (1, 1, Nz) and will contain the
        % moment1 of cloud or theta_e.
        %
        % For w-w-w, MOMENTS will be (3, 3, Nz) and the first row will
        % consist of (moment1, moment2, moment3).
        %
        % For flux, MOMENTS will be (2, 2, Nz) and the first column will
        % consist of (moment1 of the first variable, moment1 of the second
        % variable), and the first row in the second column will consist
        % of the covariance between the two variables.
        if (strcmp(Vname, 'cloud') || strcmp(Vname, 'theta_e'))
          X = 1;
          OUT_MOMENTS = squeeze(MOMENTS(1,1,:));
        elseif (strcmp(Vname, 'w-w-w'))
          OUT_MOMENTS = zeros([ 3 1 Nz ]);
          X = [ 1 2 3 ];
          OUT_MOMENTS(1,1,:) = squeeze(MOMENTS(1,1,:));
          OUT_MOMENTS(2,1,:) = squeeze(MOMENTS(1,2,:));
          OUT_MOMENTS(3,1,:) = squeeze(MOMENTS(1,3,:));
        else
          OUT_MOMENTS = zeros([ 3 1 Nz ]);
          X = [ 1 2 3 ];
          OUT_MOMENTS(1,1,:) = squeeze(MOMENTS(1,1,:));
          OUT_MOMENTS(2,1,:) = squeeze(MOMENTS(2,1,:));
          OUT_MOMENTS(3,1,:) = squeeze(MOMENTS(1,2,:));
        end

        % GenMoments() will fill levels with zero count (NPTS == 0) with nans. For
        % cloud mass we want these moments to be zero.
        if (strcmp(Vname, 'cloud'))
          OUT_MOMENTS(isnan(OUT_MOMENTS)) = 0;
        end
        
        % Write out the moment data. Put in dummy x, y and t coordinates.
        Y = 1;
        T = 1;
        
        Nx = length(X);
        OutVar = reshape(OUT_MOMENTS, [ Nx 1 Nz 1 ]);
        
        OutFile = sprintf('%s/gmd_%s_T%d_T%d_%s.h5', Ddir, Fprefix, Tstart, Tend, Case);
        fprintf('Writing: %s\n', OutFile);
        
        hdf5write(OutFile, Vname, OutVar);

        OutVar = reshape(OUT_NPTS, [ 1 1 Nz 1 ]);
        hdf5write(OutFile, 'NumPoints', OutVar, 'WriteMode', 'append');
 
        hdf5write(OutFile, 'x_coords', X, 'WriteMode', 'append');
        hdf5write(OutFile, 'y_coords', Y, 'WriteMode', 'append');
        hdf5write(OutFile, 'z_coords', Z, 'WriteMode', 'append');
        hdf5write(OutFile, 't_coords', T, 'WriteMode', 'append');
        fprintf('\n');
      end
    end
end
