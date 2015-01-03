function [ ] = GenTsPcprate(ConfigFile)

  [ Config ] = ReadConfig(ConfigFile);

  Ddir = Config.DiagDir;
  Tdir = Config. TsavgDir;

  % make sure output directory exists
  if (exist(Ddir, 'dir') ~= 7)
    mkdir(Ddir);
  end

  for icase = 1:length(Config.Cases)
    Case = Config.Cases(icase).Cname;

    fprintf('*****************************************************************\n');
    fprintf('Generating time series of pcprate:\n');
    fprintf('  Case: %s\n', Case);
    fprintf('\n');

    % Read in pcprate histograms
    InFile = sprintf('%s/hist_pcprate_%s.h5', Tdir, Case);
    InVar = '/hist_pcprate';

    HdaFile = sprintf('%s/hda_pcprate_1_%s.h5', Tdir, Case);
    HdaVar = '/hda_pcprate';

    fprintf('  Reading: %s (%s)\n', InFile, InVar);
    PHIST = squeeze(h5read(InFile, InVar));
    BINS  = squeeze(h5read(InFile, '/x_coords')); % mm/hr
    T     = squeeze(h5read(InFile, '/t_coords')) ./ 3600; % hr

    Nb = length(BINS);
    Nt = length(T);

    % Read in pcprate average
    % data is (2,t)
    %   (1,:) are sums
    %   (2,:) are counts
    fprintf('  Reading: %s (%s)\n', HdaFile, HdaVar);
    HDA_DATA = squeeze(h5read(HdaFile, HdaVar));
    HDA_PRATE_1 = HDA_DATA(1,:) ./ HDA_DATA(2,:);

    % MWIND will be: (b,t)
    % For each time step, reduce the histogram to an average
    % Try different filterings of the precip rate values
    Ptile = 50;
    PRATE     = ReduceHists(PHIST, 1, BINS, 'wtmean', Ptile);
    PRATE_0p1 = ReduceHists(PHIST(21:end,:), 1, BINS(21:end), 'wtmean', Ptile); % select rate >= 0.1 mm/h
    PRATE_1   = ReduceHists(PHIST(31:end,:), 1, BINS(31:end), 'wtmean', Ptile); % select rate >= 1.0 mm/h
    PRATE_10  = ReduceHists(PHIST(41:end,:), 1, BINS(41:end), 'wtmean', Ptile); % select rate >= 10.0 mm/h

    % Write out time series of average precip rate
    OutFile = sprintf('%s/ts_avg_pcprate_%s.h5', Ddir, Case);
    fprintf('    Writing: %s\n', OutFile)
    fprintf('\n');

    % If the file exists, remove it so that the HDF5 commands
    % can create a new file.
    if (exist(OutFile, 'file') == 2)
      delete(OutFile);
    end
    h5create(OutFile, '/avg_pcprate',  size(PRATE));
    h5create(OutFile, '/avg_pcprate_0p1',  size(PRATE_0p1));
    h5create(OutFile, '/avg_pcprate_1',  size(PRATE_1));
    h5create(OutFile, '/avg_pcprate_10',  size(PRATE_10));
    h5create(OutFile, '/hda_pcprate_1',  size(HDA_PRATE_1));
    h5create(OutFile, '/time', size(T));

    h5write(OutFile, '/avg_pcprate',  PRATE);
    h5write(OutFile, '/avg_pcprate_0p1',  PRATE_0p1);
    h5write(OutFile, '/avg_pcprate_1',  PRATE_1);
    h5write(OutFile, '/avg_pcprate_10',  PRATE_10);
    h5write(OutFile, '/hda_pcprate_1',  HDA_PRATE_1);
    h5write(OutFile, '/time', T);
  end
end 
