function [ ] = GenHistMeas( ConfigFile )
% GenHistMeas generate various measurments from histogram data
%

% Read in the config data
[ Config ] = ReadConfig(ConfigFile);

DiagDir = Config.DiagDir;

% Make sure output directory exists
if (exist(DiagDir, 'dir') ~= 7)
    mkdir(DiagDir);
end

for icase = 1:length(Config.Cases)
  Case = Config.Cases(icase).Cname;
  for ihmeas = 1: length(Config.Hmeas)
    Name = Config.Hmeas(ihmeas).Name;
    InDir = Config.Hmeas(ihmeas).InDir;
    Fprefix = Config.Hmeas(ihmeas).Fprefix;
    Vname = Config.Hmeas(ihmeas).Rvar;

    fprintf('***********************************************************************\n');
    fprintf('Generating Histogram Measurements:\n');
    fprintf('  Name: %s\n', Name);
    fprintf('  Variable: %s\n', Vname);
    fprintf('  Case: %s\n', Case);
    fprintf('\n');

    InFile  = sprintf('%s/%s_%s.h5', InDir, Fprefix, Case);
    OutFile = sprintf('%s/%s_%s.h5', DiagDir, Name, Case);
    Hdset   = sprintf('/%s', Vname);

    % Read in the histogram data. HIST will be organized as (r,b,z,t) where
    %    r --> radius
    %    b --> histogram bins
    %    z --> heights
    %    t --> time
    fprintf('Reading file: %s\n', InFile);
    fprintf('\n');
    R = hdf5read(InFile, '/x_coords');
    B = hdf5read(InFile, '/y_coords');
    Z = hdf5read(InFile, '/z_coords');
    T = hdf5read(InFile, '/t_coords');
    HIST = hdf5read(InFile, Hdset);
    [ Nr, Nb, Nz, Nt ] = size(HIST);

    % Want to reduce the bins dimension to a single number. Use
    % various measurements for this (weighted mean, max, center of
    % mass). These can then be used in a time series or a series
    % across the storm radius.

    % Do a weighted (by the histogram counts) average of the bin values.
    %   sum (bin(i) * counts(i)) / sum(counts(i)
    %
    % Numerator
    % Replicate the bin values (B) into an array (r,b,z,t) that matches
    % the size of HIST and has the bin values running along the bins
    % dimension so that an element wise multiply can be done with HIST
    % to form the bin(i) * counts(i) products everywhere. Then sum these
    % up along the bins dimension.
    BROW = reshape(B, [ 1 Nb ] );  % B needs to be a row vector
    BINS = repmat(BROW, [ Nr Nz Nt ]);
    BINS = reshape(BINS, [ Nr Nb Nz Nt ]);
    SUM_BC = squeeze(sum(HIST .* BINS, 2));

    % Denominator
    % Sum along the bins dimension which holds the counts.
    % It's possible to get all counts equal to zero since we can filter
    % out input data before making the counts. In the case where an entry
    % in SUM_C comes out to be zero, it means all the counts were zero
    % and for the purpose of this diagnostic we want the 0/0 operations
    % to result in a zero instead of NaN. MATLAB rightly results in a NaN
    % for 0/0, so we can force this to zero by changing all of the places
    % where SUM_C == 0 to a non-zero number (now were doing 0/x where x
    % is not zero).
    SUM_C = squeeze(sum(HIST,2));
    SUM_C(SUM_C == 0) = 1;

    HIST_WMEAN = SUM_BC ./ SUM_C;

    % Find the indices of the max values along the bins dimension (number 2)
    % Then use that result to create an array organized as (r,z,t) that
    % contains the bin value corresponding to the index where the max count
    % existed.
    [ MAXC_V, MAXC_I ] = max(HIST,[],2);
    HIST_MAX = squeeze(B(MAXC_I));

    % For every histogram (bins dimension: 2) find the index nearest where
    % the area under the curve is half the total area under the curve. Then
    % use these indices to place the corresponding bin values into the result.
    % This will track where the "center of mass" of the histogram distribution
    % lies.
    %
    % Accomplish this by subtracting off 1/2 the total sum from the cumulative
    % sum, taking the absolute value of the result, and locating the indices
    % where the minimum occurs (closest differece to zero). Then use these
    % indices to build a result (r,z,t) that contains the corresponding bin
    % values.
    CSUM = cumsum(HIST,2);
    SUM = sum(HIST,2);
    SUM = repmat(SUM, [ 1 Nb 1 1 ]);
    CSHALF = abs(CSUM - (SUM/2));
    [ HALF_V, HALF_I ] = min(CSHALF, [], 2);
    HIST_COM = squeeze(B(HALF_I));

    % Also compute the time mean of the histogram measurements
    % All variables are (r,z,t), this will convert to (r,z).
    HIST_WMEAN_TMEAN = mean(HIST_WMEAN,3);
    HIST_MAX_TMEAN   = mean(HIST_MAX,3);
    HIST_COM_TMEAN   = mean(HIST_COM,3);

    % If doing 'tempc'
    % Calculate a line representing the freezing level
    % Don't use >= 0 in the find command because the wtmean method
    % can result in multiple zeros in a profile
    if (strcmp(Vname, 'tempc'))
      FRZLEV_WMEAN = zeros(1,Nr);
      FRZLEV_MAX = zeros(1,Nr);
      FRZLEV_COM = zeros(1,Nr);
      for ir = 1:Nr
        i = find(HIST_WMEAN_TMEAN(ir,:) > 0, 1, 'last');
        if (length(i) == 0)
          i = 1;
        end
        FRZLEV_WMEAN(ir) = Z(i);
        i = find(HIST_MAX_TMEAN(ir,:) > 0, 1, 'last');
        if (length(i) == 0)
          i = 1;
        end
        FRZLEV_MAX(ir)   = Z(i);
        i = find(HIST_COM_TMEAN(ir,:) > 0, 1, 'last');
        if (length(i) == 0)
          i = 1;
        end
        FRZLEV_COM(ir)   = Z(i);
      end
    end

    fprintf('Writing file: %s\n', OutFile);
    fprintf('\n');

    Hdset = sprintf('/%s_wtmean', Vname);
    hdf5write(OutFile, Hdset, HIST_WMEAN);
    Hdset = sprintf('/%s_max', Vname);
    hdf5write(OutFile, Hdset, HIST_MAX, 'WriteMode', 'append');
    Hdset = sprintf('/%s_com', Vname);
    hdf5write(OutFile, Hdset, HIST_COM, 'WriteMode', 'append');

    Hdset = sprintf('/%s_wtmean_tmean', Vname);
    hdf5write(OutFile, Hdset, HIST_WMEAN_TMEAN, 'WriteMode', 'append');
    Hdset = sprintf('/%s_max_tmean', Vname);
    hdf5write(OutFile, Hdset, HIST_MAX_TMEAN, 'WriteMode', 'append');
    Hdset = sprintf('/%s_com_tmean', Vname);
    hdf5write(OutFile, Hdset, HIST_COM_TMEAN, 'WriteMode', 'append');

    if (strcmp(Vname, 'tempc'))
      Hdset = sprintf('/frzlev_wtmean');
      hdf5write(OutFile, Hdset, FRZLEV_WMEAN, 'WriteMode', 'append');
      Hdset = sprintf('/frzlev_max');
      hdf5write(OutFile, Hdset, FRZLEV_MAX, 'WriteMode', 'append');
      Hdset = sprintf('/frzlev_com');
      hdf5write(OutFile, Hdset, FRZLEV_COM, 'WriteMode', 'append');
    end

    hdf5write(OutFile, '/x_coords', R, 'WriteMode', 'append');
    hdf5write(OutFile, '/y_coords', B, 'WriteMode', 'append');
    hdf5write(OutFile, '/z_coords', Z, 'WriteMode', 'append');
    hdf5write(OutFile, '/t_coords', T, 'WriteMode', 'append');

    fprintf('\n');
  end
end

end
