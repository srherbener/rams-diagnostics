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

%for ismeas = 1: length(Config.Smeas)
for ismeas = 1: 1
  % clear these out because they may need to change size, and without the clear
  % they will just retain the largest size (and associated data) that they have
  % been set to so far
  clear Xvals;
  clear RDATA;

  Name = Config.Smeas(ismeas).Name;
  InDir = Config.Smeas(ismeas).InDir;
  Fprefix = Config.Smeas(ismeas).Fprefix;
  Vname = Config.Smeas(ismeas).Rvar;
  Xgroup = Config.Smeas(ismeas).Xgroup;
  Ygroup = Config.Smeas(ismeas).Ygroup;
  Xvals = Config.Smeas(ismeas).Xvals'; % note transpose to get dimensions matching Yvals below

  fprintf('***********************************************************************\n');
  fprintf('Generating Slope time series:\n');
  fprintf('  Name: %s\n', Name);
  fprintf('  Variable: %s\n', Vname);
  fprintf('  Bin group size in x: %d\n', Xgroup);
  fprintf('  Bin group size in y: %d\n', Ygroup);
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
      %    y --> LTSS bins
      %    z --> counts: 1 -> Nt, 2 -> Nr
      %    t --> time
      fprintf('  Reading file: %s, Dataset: %s\n', InFile, Hdset);
      fprintf('\n');
      X = hdf5read(InFile, '/x_coords');
      Y = hdf5read(InFile, '/y_coords');
      T = hdf5read(InFile, '/t_coords');
      COUNTS = hdf5read(InFile, Hdset);

      % Generate the bins for the slope calculations. Xgroup (and Ygroup) determines
      % how to reduce bins (by combining the counts for adjacent groups of bins).
      % If Xgroup is 5, then every 5 adjacent bins have the counts summed
      % together to form one output bin, so 100 bins would get reduced to 20 bins,
      % where the first output bin would consitst of the counts summed together
      % from the input bins 1-5, second output bins from the input bins 6-10, etc.
      [ NR, NT, XL, XU, YL, YU ] = GenSlopeBins(COUNTS, X, Y, Xgroup, Ygroup);
NR
NT

      % At this point we have:
      %    NR(x,y,t) - count of grid cells that are raining, by bin
      %    NT(x,y,t) - count of total grid cells, by bin
      %    XL(x)   - lower edge of each x bin
      %    XU(x)   - upper edge of each x bin
      %    YL(y)   - lower edge of each y bin
      %    YU(y)   - upper edge of each y bin
      RDATA(:,:,:,icase) = NR ./ NT;
    end

    % RDATA is now orgainized as (x,y,t,c) where x are the LWP bins, y are the LTSS bins,
    % t are the times and c are the cases. If there were no cells selected for a particular
    % LWP,LTSS bin, then both the Nr and Nt values were zero, and RDATA contains a nan for
    % the 0/0 operation. Want these to be zeros so change all the nan's into zeros.
    %
    % WARNING: this assumes that the counts got constructed correctly. Ie, you can't have a
    % non-zero Nr without a non-zero Nt, so the only time you get a nan is when Nr = Nt = 0.
    RDATA(isnan(RDATA)) = 0;
    [ Nx, Ny, Nt, Nc ] = size(RDATA);

    % Now go through each set of ratios for each LWP bin, and calculate the slope (regression
    % coefficient) relating the given variable (Xvals) to the ratios.
    SLOPES     = zeros([ Nx Ny Nt ]);
    YINTS      = zeros([ Nx Ny Nt ]);
    CORCOEFFS  = zeros([ Nx Ny Nt ]);
    for it = 1:Nt % loop over time
      for ix = 1:Nx % loop over LWP bins
        for iy = 1:Ny % loop over LTSS bins
          Yvals = squeeze(RDATA(ix,iy,it,:));
          [ SLOPES(ix,iy,it) YINTS(ix,iy,it) CORCOEFFS(ix,iy,it) ] = RegFit(Xvals, Yvals); 
        end
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

    hdf5write(OutFile, '/LwpBinLower', XL, 'WriteMode', 'append');
    hdf5write(OutFile, '/LwpBinUpper', XU, 'WriteMode', 'append');
    hdf5write(OutFile, '/LtssBinLower', YL, 'WriteMode', 'append');
    hdf5write(OutFile, '/LtssBinUpper', YU, 'WriteMode', 'append');
    hdf5write(OutFile, '/Time', T, 'WriteMode', 'append');

    fprintf('\n');
  end
end

end
