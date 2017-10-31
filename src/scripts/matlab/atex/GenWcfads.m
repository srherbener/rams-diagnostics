function [ ] = GenWcfads(ConfigFile)
% GenWcfads generate cfads of w data averaged over a time interval

    % Read the config file to get the structure of how the data is laid out in
    % the file system.
    [ Config ] = ReadConfig(ConfigFile);

    Tdir = Config.TsavgDir;
    Ddir = Config.DiagDir;

    Wname = 'hist_w';

    Wfcases = {
     'hist_w'
%     'hist_w_nz'
     };
    Ofcases = {
     'cfad_w'
%     'cfad_w_nz'
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

    for icase = 1:length(Config.Cases)
        Case = Config.Cases(icase).Cname;

        for ifc = 1:length(Wfcases)
            Wfcase = Wfcases{ifc};
            Ofcase = Ofcases{ifc};

            Wfile   = sprintf('%s/%s_%s.h5', Tdir, Wfcase, Case);
    
            fprintf('***************************************************************\n');
            fprintf('Generating w cfad:\n');
            fprintf('  Case: %s\n', Case);
            fprintf('  Input w histogram file: %s\n', Wfile);
            fprintf('    Var name: %s\n', Wname);
            fprintf('\n');

            % Read in data
            %  NOTE: T is in seconds and Tstart,Tend are in hours
            W_HIST = squeeze(hdf5read(Wfile, Wname));
            X = squeeze(hdf5read(Wfile, '/x_coords'));
            Y = 1; % Y is a dummy variable
            Z = squeeze(hdf5read(Wfile, '/z_coords'));
            T = squeeze(hdf5read(Wfile, '/t_coords')) / 3600;  % hours

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

              OUT_W_CFAD     = reshape(W_CFAD,     [ Nx 1 Nz 1 ]);
              fprintf('  Writing: %s\n', OutFile);
              hdf5write(OutFile, '/W_CFAD',     OUT_W_CFAD);
              hdf5write(OutFile, '/x_coords',   X,           'WriteMode', 'append');
              hdf5write(OutFile, '/y_coords',   Y,           'WriteMode', 'append');
              hdf5write(OutFile, '/z_coords',   Z,           'WriteMode', 'append');
              hdf5write(OutFile, '/t_coords',   Tdummy,      'WriteMode', 'append');
              fprintf('\n');
            end
        end
    end
end
