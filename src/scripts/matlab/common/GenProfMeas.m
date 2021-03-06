function [ ] = GenProfMeas( ConfigFile )
% GenProfMeas generate profile data: CFAD, Time series of representative profile
%
% CFAD
%   This routine will select histogram data from given radius range
%   and time range. Then the histogram count data will simply be
%   summed up over the selected space. This will result in 2D
%   data of the form: (b,z) where b is the histogram bins and
%   z is height.
%
% Time Series of Representative Profile
%   This routine will select histogram data from the give radius range
%   only. Then the histogram data will be summed up over the selected
%   space, and then reduced using the given method. This will result
%   in a time series of a vertical profile.
%

% Read in the config data
[ Config ] = ReadConfig(ConfigFile);

DiagDir = Config.DiagDir;
UndefVal = Config.UndefVal;

% Make sure output directory exists
if (exist(DiagDir, 'dir') ~= 7)
    mkdir(DiagDir);
end

for icase = 1:length(Config.Cases)
  Case = Config.Cases(icase).Cname;
  for icfad = 1: length(Config.Pmeas)
    Name = Config.Pmeas(icfad).Name;
    InDir = Config.Pmeas(icfad).InDir;
    Fprefix = Config.Pmeas(icfad).Fprefix;
    Vname = Config.Pmeas(icfad).Rvar;
    Method = Config.Pmeas(icfad).Method;
    Rmin = Config.Pmeas(icfad).Rmin;
    Rmax = Config.Pmeas(icfad).Rmax;
    Tmin = Config.Pmeas(icfad).Tmin;
    Tmax = Config.Pmeas(icfad).Tmax;

    fprintf('***********************************************************************\n');
    fprintf('Generating CFAD and profile time series:\n');
    fprintf('  Name: %s\n', Name);
    fprintf('  Variable: %s\n', Vname);
    fprintf('  Method: %s\n', Method);
    fprintf('  Case: %s\n', Case);
    fprintf('  Rmin, Rmax: %.2f, %.2f\n', Rmin, Rmax);
    fprintf('  Tmin, Tmax: %.2f, %.2f\n', Tmin, Tmax);
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

    % Select the space to do the summing
    R1 = find(R/1000 >= Rmin, 1, 'first'); % Rmin, Rmax are in km
    R2 = find(R/1000 <= Rmax, 1, 'last');
    T1 = find(T/3600 >= Tmin, 1, 'first'); % Tmin and Tmax are in hrs
    T2 = find(T/3600 <= Tmax, 1, 'last');

    if (strcmp(Method, 'azavg'))
      % Will be workign with azimuthally averaged data, so there are no histograms for the 
      % horizontal data. Can make PROF, PROF_TS and PROF_RS but cannot make CFAD. Just use
      % PROF for CFAD so that the output section has something to write out.
      %
      % HIST will be organized as (r,d,z,t) where d ia a dummy dimension (size 1, value 0)

      % In cases where no data was selected, the resulting entry will be UndefVal. Change
      % these to nans before doing averaging.
      HIST(HIST == UndefVal) = nan;
      
      % Single profile
      PROF = squeeze(nanmean(nanmean(HIST(R1:R2,:,:,T1:T2),1),4));
      CFAD = PROF;

      % Profile time series (take mean over radius)
      PROF_TS = squeeze(nanmean(HIST(R1:R2,:,:,T1:T2),1));

      % Profile radial series (take mean over time)
      PROF_RS = squeeze(nanmean(HIST(R1:R2,:,:,T1:T2),4));
    else
      % Will be working on histogram data, so create profiles with ReduceHists()

      % CFAD
      % Sum across the r and t dimensions (1st and 4th)
      % limiting to the selected values
      CFAD = squeeze(sum(sum(HIST(R1:R2,:,:,T1:T2),1),4));
  
      % Single Profile (corresponiding to CFAD)
      [ PROF ] = ReduceHists(CFAD, 1, B, Method);
  
      % Profile Time Series
      %
      % Sum across only the r dimension, then reduce the bin dimension.
      R_HISTS = squeeze(sum(HIST(R1:R2,:,:,T1:T2),1));
      [ PROF_TS ] = ReduceHists(R_HISTS, 1, B, Method);
  
      % Profile Radial Series
      %
      % Sum acroos only the t dimension, then reduce the bin dimension.
      T_HISTS = squeeze(sum(HIST(R1:R2,:,:,T1:T2),4));
      [ PROF_RS ] = ReduceHists(T_HISTS, 2, B, Method);
    end

    fprintf('Writing file: %s\n', OutFile);
    fprintf('\n');

    Hdset = sprintf('/Cfad_%s', Vname);
    hdf5write(OutFile, Hdset, CFAD);
    Hdset = sprintf('/Prof_%s', Vname);
    hdf5write(OutFile, Hdset, PROF, 'WriteMode', 'append');
    Hdset = sprintf('/ProfTs_%s', Vname);
    hdf5write(OutFile, Hdset, PROF_TS, 'WriteMode', 'append');
    Hdset = sprintf('/ProfRs_%s', Vname);
    hdf5write(OutFile, Hdset, PROF_RS, 'WriteMode', 'append');

    hdf5write(OutFile, '/Radius', R(R1:R2), 'WriteMode', 'append');
    hdf5write(OutFile, '/Bins', B, 'WriteMode', 'append');
    hdf5write(OutFile, '/Height', Z, 'WriteMode', 'append');
    hdf5write(OutFile, '/Time', T(T1:T2), 'WriteMode', 'append');

    fprintf('\n');
  end
end

end
