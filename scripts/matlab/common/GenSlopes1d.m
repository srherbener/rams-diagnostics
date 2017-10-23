function [ ] = GenSlopes1d(ConfigFile)
% GenSlopes1d generate slope data from POP diagnostic

% Read the config file to get the structure of how the data is laid out in
% the file system.
[ Config ] = ReadConfig(ConfigFile);

Ddir = Config.DiagDir;

% Make sure output directory exists
if (exist(Ddir, 'dir') ~= 7)
    mkdir(Ddir);
end

for ismeas = 1: length(Config.Smeas1d)
  % clear these out because they may need to change size, and without the clear
  % they will just retain the largest size (and associated data) that they have
  % been set to so far
  clear Xvals;
  clear RDATA;

  Name = Config.Smeas1d(ismeas).Name;
  InDir = Config.Smeas1d(ismeas).InDir;
  Fprefix = Config.Smeas1d(ismeas).Fprefix;
  Vname = Config.Smeas1d(ismeas).Rvar;
  Xmin = Config.Smeas1d(ismeas).Xmin;
  Xmax = Config.Smeas1d(ismeas).Xmax;
  Xgroup = Config.Smeas1d(ismeas).Xgroup;
  Xvals = Config.Smeas1d(ismeas).Xvals'; % note transpose to get dimensions matching Yvals below

  fprintf('***********************************************************************\n');
  fprintf('Generating Slope time series:\n');
  fprintf('  Name: %s\n', Name);
  fprintf('  Variable: %s\n', Vname);
  fprintf('  X bin selection:\n');
  fprintf('    Min: %.2f\n', Xmin);
  fprintf('    Max: %.2f\n', Xmax);
  fprintf('    Bin group size: %d\n', Xgroup);
  fprintf('  X values for the regression fit:\n');
  for ix = 1:length(Xvals)
    fprintf('    Xvals(%d): %f\n', ix, Xvals(ix));
  end
  fprintf('\n');

  OutFile = sprintf('%s/%s.h5', Ddir, Name);

  ips = Config.Smeas1d(ismeas).PSnum;
  if (ips == 0)
    fprintf('  WARNING: skipping Smeas1d number %d due to no associated PlotSet\n', ismeas)
  else
    % For each case in the plot set, read in the count data, form the ratio between
    % the two sets of counts (y(2) / y(1)).
    for icase = 1:Config.PlotSets(ips).Ncases
      Case = Config.PlotSets(ips).Cases(icase).Cname;

      InFile  = sprintf('%s/%s_%s.h5', InDir, Fprefix, Case);
      Hdset   = sprintf('/%s', Vname);

      % Read in the count data. COUNTS will be organized as (x,y,z,t) where
      %    x --> LWP bins
      %    y --> LTSS bins
      %    z --> counts: 1 -> Nt, 2 -> Nr
      %    t --> time
      fprintf('  Reading file: %s, Dataset: %s\n', InFile, Hdset);
      fprintf('\n');
      X = hdf5read(InFile, '/x_coords');
      T = hdf5read(InFile, '/t_coords');
      COUNTS = hdf5read(InFile, Hdset);
      
      % The y-dimension (LTSS) was set up to place everything into one bin
      % so sum up across those bins. COUNTS will then be (x,z,t).
      COUNTS = squeeze(sum(COUNTS,2));

      % Generate the bins for the slope calculations. Xgroup determines
      % how to reduce bins (by combining the counts for adjacent groups of bins).
      % If Xgroup is 5, then every 5 adjacent bins have the counts summed
      % together to form one output bin, so 100 bins would get reduced to 20 bins,
      % where the first output bin would consitst of the counts summed together
      % from the input bins 1-5, second output bins from the input bins
      % 6-10, etc.
      [ NR, XL, XU ] = GenCountBins1d(COUNTS(:,2,:), X, Xmin, Xmax, Xgroup);
      [ NT, XL, XU ] = GenCountBins1d(COUNTS(:,1,:), X, Xmin, Xmax, Xgroup);

      % At this point we have:
      %    NR(x,t) - count of grid cells that are raining, by bin
      %    NT(x,t) - count of total grid cells, by bin
      %    XL(x)   - lower edge of each x bin
      %    XU(x)   - upper edge of each x bin
      RDATA(:,:,icase) = NR ./ NT;
    end

    % RDATA is now orgainized as (x,t,c) where x are the LWP bins,
    % t are the times and c are the cases. If there were no cells selected for a particular
    % LWP bin, then both the Nr and Nt values were zero, and RDATA contains a nan for
    % the 0/0 operation. Want these to be zeros so change all the nan's into zeros.
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
    NUSED      = zeros([ Nx Nt ]);
    for it = 1:Nt % loop over time
      for ix = 1:Nx % loop over LWP bins
        Yvals = squeeze(RDATA(ix,it,:));
        [ SLOPES(ix,it) YINTS(ix,it) CORCOEFFS(ix,it) NUSED(ix,it) ] = RegFit(Xvals, Yvals); 
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
    Hdset = sprintf('/Nused_%s', Vname);
    hdf5write(OutFile, Hdset, NUSED, 'WriteMode', 'append');

    hdf5write(OutFile, '/LwpBinLower', XL, 'WriteMode', 'append');
    hdf5write(OutFile, '/LwpBinUpper', XU, 'WriteMode', 'append');
    hdf5write(OutFile, '/Time', T, 'WriteMode', 'append');

    fprintf('\n');
  end
end

end
