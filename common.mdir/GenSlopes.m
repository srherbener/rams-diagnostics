function [ ] = GenSlopes(ConfigFile)
% GenSlopes generate slope data from POP diagnostic

% Read the config file to get the structure of how the data is laid out in
% the file system.
[ Config ] = ReadConfig(ConfigFile);

Tdir = Config.TsavgDir;
Ddir = Config.DiagDir;

% Make sure output directory exists
if (exist(Ddir, 'dir') ~= 7)
    mkdir(Ddir);
end

for ismeas = 1: length(Config.Smeas)
  % clear these out because they may need to change size, and without the clear
  % they will just retain the largest size (and associated data) that they have
  % been set to so far
  clear Xvals;
  clear RDATA;

  Name = Config.Smeas(ismeas).Name;
  InDir = Config.Smeas(ismeas).InDir;
  Fprefix = Config.Smeas(ismeas).Fprefix;
  Vname = Config.Smeas(ismeas).Rvar;
  Xvals = Config.Smeas(ismeas).Xvals'; % note transpose to get dimensions matching Yvals below

  fprintf('***********************************************************************\n');
  fprintf('Generating Slope time series:\n');
  fprintf('  Name: %s\n', Name);
  fprintf('  Variable: %s\n', Vname);
  fprintf('  X values for the regression fit:\n');
  for ix = 1:length(Xvals)
    fprintf('    Xvals(%d): %f\n', ix, Xvals(ix));
  end
  fprintf('\n');

  OutFile = sprintf('%s/%s.h5', Ddir, Name);

  ips = Config.Smeas(ismeas).PSnum;
  if (ips == 0)
    fprintf('  WARNING: skipping Smeas number %d due to no associated PlotSet\n', ismeas)
  else
    % For each case in the plot set, read in the count data, form the ratio between
    % the two sets of counts (y(2) / y(1)).
    for icase = 1:Config.PlotSets(ips).Ncases
      Case = Config.PlotSets(ips).Cases(icase).Cname;

      InFile  = sprintf('%s/%s_%s.h5', InDir, Fprefix, Case);
      Hdset   = sprintf('/%s', Vname);

      % Read in the count data. COUNTS will be organized as (x,y,z,t) where
      %    x --> LWP bins
      %    y --> counts: 1: Nt, 2: Nr
      %    z --> dummy
      %    t --> time
      fprintf('  Reading file: %s, Dataset: %s\n', InFile, Hdset);
      fprintf('\n');
      X = hdf5read(InFile, '/x_coords');
      T = hdf5read(InFile, '/t_coords');
      COUNTS = squeeze(hdf5read(InFile, Hdset));

      % Do a bit-wise divide of the 2st column of y by the 1st column of y. Ie, form
      % the values Nr/Nt. This will reduce y to size-1 dimension so squeeze it out -->
      % results in the array being organized as (x,t) where x are the Nr/Nt ratios.
      % Keep a copy of the ratios for each case.
      RDATA(:,:,icase) = squeeze(COUNTS(:,2,:) ./ COUNTS(:,1,:));
    end

    % RDATA is now orgainized as (x,t,c) where x are the LWP bins, t are the times and
    % c are the cases. If there were no cells selected for a particular LWP bin, then both
    % the Nr and Nt values were zero, and RDATA contains a nan for the 0/0 operation. Want
    % these to be zeros so change all the nan's into zeros.
    %
    % WARNING: this assumes that the counts got constructed correctly. Ie, you can't have a
    % non-zero Nr without a non-zero Nt, so the only time you get a nan is when Nr = Nt = 0.
    RDATA(isnan(RDATA)) = 0;
    [ Nx, Nt, Nc ] = size(RDATA);

    % Now go through each set of ratios for each LWP bin, and calculate the slope (regression
    % coefficient) relating the given variable (Xvals) to the ratios.
    SLOPES     = zeros([ Nx Nt ]);
    YINTS      = zeros([ Nx Nt ]);
    CORCOEFFS  = zeros([ Nx Nt ]);
    for it = 1:size(RDATA, 2) % loop over time
      for ix = 1:size(RDATA, 1) % loop over LWP bins
        Yvals = squeeze(RDATA(ix,it,:));
        [ SLOPES(ix,it) YINTS(ix,it) CORCOEFFS(ix,it) ] = RegFit(Xvals, Yvals); 
      end
    end

    fprintf('Writing file: %s\n', OutFile);
    fprintf('\n');

    Hdset = sprintf('/Slopes_%s', Vname);
    hdf5write(OutFile, Hdset, SLOPES);
    Hdset = sprintf('/Yints_%s', Vname);
    hdf5write(OutFile, Hdset, YINTS, 'WriteMode', 'append');
    Hdset = sprintf('/CorCoeffs_%s', Vname);
    hdf5write(OutFile, Hdset, CORCOEFFS, 'WriteMode', 'append');

    hdf5write(OutFile, '/LwpBins', X, 'WriteMode', 'append');
    hdf5write(OutFile, '/Time', T, 'WriteMode', 'append');

    fprintf('\n');
  end
end

end
