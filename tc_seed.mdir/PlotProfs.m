function [ ] = PlotProfs(ConfigFile)
% PlotProfs generate vertical profile plots

% Read the config file to get the structure of how the data is laid out in
% the file system.
[ Config ] = ReadConfig(ConfigFile);

Ddir = Config.DiagDir;
Pdir = Config.PlotDir;

% Find and replace underscores in Ptitle, Ylabel with blank spaces
for iplot = 1:length(Config.ProfPlots)
    clear Profs;
    clear LegText;
    clear Lspecs;

    Fprefix = Config.ProfPlots(iplot).Fprefix;
    Var = Config.ProfPlots(iplot).Var;

    % config for axes
    Xstart = Config.ProfPlots(iplot).Xspec(1);
    Xinc   = Config.ProfPlots(iplot).Xspec(2);
    Xend   = Config.ProfPlots(iplot).Xspec(3);
    Xvals  = (Xstart:Xinc:Xend);

    Zmin = Config.ProfPlots(iplot).Zmin;
    Zmax = Config.ProfPlots(iplot).Zmax;
    
    Ptype = Config.ProfPlots(iplot).Type;

    % If doing a diff plot, read in the control profile
    if (~isempty(strfind(Ptype, 'diff')))
      Case = Config.ControlCase;
      Hfile = sprintf('%s/%s_%s.h5', Ddir, Fprefix, Case);
      Hdset = sprintf('/Prof_%s', Var);
      fprintf('Reading Control Case: %s\n', Case);
      fprintf('  HDF5 file: %s, dataset: %s\n', Hfile, Hdset);
      PROF_CNTL = squeeze(hdf5read(Hfile, Hdset));
    end

    Ptitle = sprintf('%s', Config.ProfPlots(iplot).Title);
    Xlabel = Config.ProfPlots(iplot).Xlabel;
    Zlabel = Config.ProfPlots(iplot).Zlabel;
    LegLoc = Config.ProfPlots(iplot).LegLoc;
    OutFile = sprintf('%s/%s', Pdir, Config.ProfPlots(iplot).OutFile);
    
    % make sure output directory exists
    if (exist(Pdir, 'dir') ~= 7)
        mkdir(Pdir);
    end

    iprof = 0;
    ips = Config.ProfPlots(iplot).PSnum;
    if (ips == 0)
      fprintf('WARNING: skipping ProfPlot number %d due to no associated PlotSet\n', iplot)
    else
      for icase = 1:Config.PlotSets(ips).Ncases
        Case = Config.PlotSets(ips).Cases(icase).Cname;
        if (strcmp(Case, Config.ControlCase) && (~isempty(strfind(Ptype, 'diff'))))
          % skip the control case when doing a diff plot
          continue; 
        end

        iprof = iprof + 1;
        LegText(iprof) = { Config.PlotSets(ips).Cases(icase).Legend };
        Lspecs(iprof) = { Config.PlotSets(ips).Cases(icase).Lspec };

        % Var is organized (b,z) in the file (CFAD).
        % Bins holds the bin values, Height holds the z values.
        Hfile = sprintf('%s/%s_%s.h5', Ddir, Fprefix, Case);
        Hdset = sprintf('/Prof_%s', Var);
        fprintf('  HDF5 file: %s, dataset: %s\n', Hfile, Hdset);
        PROF_EXP = squeeze(hdf5read(Hfile, Hdset));
        Z = hdf5read(Hfile, 'Height')/1000; % km

        % Trim off the selected z range from the lhv data
        % Each profile goes into a row of LHV
        Z1 = find(Z >= Zmin, 1, 'first');
        Z2 = find(Z <= Zmax, 1, 'last');
        Zvals = Z(Z1:Z2);
        if (~isempty(strfind(Ptype, 'diff')))
          DIFF = PROF_EXP(Z1:Z2) - PROF_CNTL(Z1:Z2);
          
          if (strcmp(Ptype, 'pdiff'))
              DIFF = (DIFF ./ abs(PROF_CNTL(Z1:Z2))) * 100;
          end
          
          Profs(iprof,:) = DIFF;
        else
          Profs(iprof,:) = PROF_EXP(Z1:Z2);
        end
      end
    end

    fprintf('Writing plot file: %s\n', OutFile);
    PlotProfSet(Xvals, Zvals, Profs, Xlabel, Zlabel, Ptitle, Lspecs, LegText, LegLoc, OutFile);
    fprintf('\n');
end
