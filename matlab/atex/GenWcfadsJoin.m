function [ ] = GenWcfadsJoin(ConfigFile)
% GenWcfadsJoin generate cfads of w data over a time interval, join two input files together to create the cfad

    % Read the config file to get the structure of how the data is laid out in
    % the file system.
    [ Config ] = ReadConfig(ConfigFile);

    Tdir = Config.TsavgDir;
    Ddir = Config.DiagDir;

    Wname = 'hist_w';

    UPfcases = {
     'hist_up_0p01'
     'hist_up_0p10'
     };
    DNfcases = {
     'hist_dn_0p01'
     'hist_dn_0p10'
     };
    Ofcases = {
     'cfad_w_0p01'
     'cfad_w_0p10'
     };

    TstartTimes = [
     12
     24
     35
     ];
    TendTimes = [
     13
     25
     36
     ];

    GroupSize = 20;

    for icase = 1:length(Config.Cases)
        Case = Config.Cases(icase).Cname;

        for ifc = 1:length(UPfcases)
            UPfcase = UPfcases{ifc};
            DNfcase = DNfcases{ifc};
            Ofcase = Ofcases{ifc};

            UPfile   = sprintf('%s/%s_%s.h5', Tdir, UPfcase, Case);
            DNfile   = sprintf('%s/%s_%s.h5', Tdir, DNfcase, Case);
    
            fprintf('***************************************************************\n');
            fprintf('Generating w cfad:\n');
            fprintf('  Case: %s\n', Case);
            fprintf('  Input updraft histogram file: %s\n', UPfile);
            fprintf('    Var name: %s\n', Wname);
            fprintf('  Input downdraft histogram file: %s\n', DNfile);
            fprintf('    Var name: %s\n', Wname);
            fprintf('\n');

            % Read in data
            %  NOTE: T is in seconds and Tstart,Tend are in hours
            UP_HIST = squeeze(hdf5read(UPfile, Wname));
            UP_X = squeeze(hdf5read(UPfile, '/x_coords'));
            DN_HIST = squeeze(hdf5read(DNfile, Wname));
            DN_X = squeeze(hdf5read(DNfile, '/x_coords'));
            Z = squeeze(hdf5read(DNfile, '/z_coords'));
            T = squeeze(hdf5read(DNfile, '/t_coords')) / 3600;  % hours

            % Combine histogram bins in the updraft and downdraft data, then join the two to create on histogram.
            % The last bin of the downdraft histogram is for w = 0. This bin will always contain a count of zero
            % due to the filtering out of zero valued vertical velocity. Therefore, throw away this bin.
            fprintf('  Combining histogram bins, group size: %d\n', GroupSize);
            fprintf('\n');

            UP = CombineBins(UP_HIST, 1, GroupSize);
            DN = CombineBins(DN_HIST(1:end-1,:,:), 1, GroupSize); % strip off the last bin (w == zero)
            W_HIST = cat(1, DN, UP);
            DN_X_TEMP = DN_X(1:end-1);  % strip off the last bin (w == zero)
            X = cat(1, DN_X_TEMP(1:GroupSize:end), UP_X(1:GroupSize:end));

            for it = 1:length(TstartTimes)
              Tstart = TstartTimes(it);
              Tend   = TendTimes(it);
              OutFile = sprintf('%s/%s_T%d_T%d_%s.h5', Ddir, Ofcase, Tstart, Tend, Case);

              fprintf('  Time Average Interval:\n');
              fprintf('    Tstart: %d (h)\n', Tstart);
              fprintf('    Tend: %d (h)\n', Tend);
              fprintf('\n');

              T1 = find(T >= Tstart, 1, 'first');
              T2 = find(T <= Tend,   1, 'last');

              % sum up histogram counts over the interval T1 to T2
              % make each level a PDF of W
              W_CFAD = squeeze(sum(W_HIST(:,:,T1:T2), 3));
              [ Nx Nz ] = size(W_CFAD);
              S = sum(W_CFAD, 1);
              S = repmat(S, [ Nx 1 ]);
              W_CFAD = W_CFAD ./ S;

              % output --> Use REVU format, 4D var, *_coords
              Tdummy = 1;  % due to summing over time, t becomes a dummy variable
              Ydummy = 1; % Y is a dummy variable

              OUT_W_CFAD     = reshape(W_CFAD,     [ Nx 1 Nz 1 ]);
              fprintf('  Writing: %s\n', OutFile);
              hdf5write(OutFile, '/W_CFAD',     OUT_W_CFAD);
              hdf5write(OutFile, '/x_coords',   X,           'WriteMode', 'append');
              hdf5write(OutFile, '/y_coords',   Ydummy,      'WriteMode', 'append');
              hdf5write(OutFile, '/z_coords',   Z,           'WriteMode', 'append');
              hdf5write(OutFile, '/t_coords',   Tdummy,      'WriteMode', 'append');
              fprintf('\n');
            end
        end
    end
end
