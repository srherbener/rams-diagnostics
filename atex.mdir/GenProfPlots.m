function [ ] = GenProfPlots(ConfigFile)
% GenProfPlots generate vertical profile plots

% Read the config file to get the structure of how the data is laid out in
% the file system.
[ Config ] = ReadConfig(ConfigFile);
    
Pname   = Config.ExpName;
UndefVal = Config.UndefVal;

Tdir = Config.TsavgDir;
Pdir = Config.PlotDir;

% Find and replace underscores in Ptitle, Ylabel with blank spaces
for iplot = 1:length(Config.ProfPlots)
    clear LHV;
    clear LegText;
    clear Lstyles;
    clear Lgscales;

    Fprefix = Config.ProfPlots(iplot).Fprefix;
    Var = Config.ProfPlots(iplot).Var;

    % config for axes
    Xstart = Config.ProfPlots(iplot).Xspec(1);
    Xinc   = Config.ProfPlots(iplot).Xspec(2);
    Xend   = Config.ProfPlots(iplot).Xspec(3);
    X = (Xstart:Xinc:Xend);

    Zmin = Config.ProfPlots(iplot).Zmin;
    Zmax = Config.ProfPlots(iplot).Zmax;

    % If doing a diff plot, read in the control profile
    Ptype = Config.ProfPlots(iplot).Type;
    if (strcmp(Ptype, 'diff'))
      Case = Config.ControlCase;
      Hfile = sprintf('%s/%s_%s.h5', Tdir, Fprefix, Case);
      fprintf('Reading Control Case: %s\n', Case);
      fprintf('  HDF5 file: %s\n', Hfile);
      LHV_DOMAVG = squeeze(hdf5read(Hfile, Var));
      LHV_DOMAVG(LHV_DOMAVG == UndefVal) = nan;
      LHV_CONTROL_ALLZ = nanmean(LHV_DOMAVG,2); % time average
    end

    Ptitle   = Config.ProfPlots(iplot).Title.Main;
    Pmarkers = Config.ProfPlots(iplot).Title.Pmarkers;
    Xlabel   = Config.ProfPlots(iplot).Xlabel;
    Zlabel   = Config.ProfPlots(iplot).Zlabel;
    LegLoc   = Config.ProfPlots(iplot).LegLoc;
    OutFile  = sprintf('%s/%s', Pdir, Config.ProfPlots(iplot).OutFile);
    
    % make sure output directory exists
    if (exist(Pdir, 'dir') ~= 7)
        mkdir(Pdir);
    end

    ihist = 0;
    ips = Config.ProfPlots(iplot).PSnum;
    if (ips == 0)
      fprintf('WARNING: skipping ProfPlot number %d due to no associated PlotSet\n', iplot)
    else
      for icase = 1:Config.PlotSets(ips).Ncases
        Case = Config.PlotSets(ips).Cases(icase).Cname;
        LegText(icase) = { Config.PlotSets(ips).Cases(icase).Legend };
        Lstyles(icase) = { Config.PlotSets(ips).Cases(icase).Lstyle };
        Lgscales(icase) = Config.PlotSets(ips).Cases(icase).Lgscale;

        % Var is organized (x,y,z,t) in the file, however x and y
        % dimension sizes are both 1. After running squeeze(), Var
        % will be reduced to (z,t).
        Hfile = sprintf('%s/%s_%s.h5', Tdir, Fprefix, Case);
        fprintf('Reading HDF5 file: %s\n', Hfile);
        LHV_DOMAVG = squeeze(hdf5read(Hfile, Var));
        LHV_DOMAVG(LHV_DOMAVG == UndefVal) = nan;
        LHV_ALLZ = nanmean(LHV_DOMAVG,2); % time average

        % Grab the heights on the first case. The z coordinate
        % values should be the same for all cases.
        if (icase == 1)
          Zall = hdf5read(Hfile, 'z_coords');
          iz1 = find(Zall >= Zmin, 1, 'first');
          iz2 = find(Zall <= Zmax, 1, 'last');
          Z = Zall(iz1:iz2);
        end

        % Trim off the selected z range from the lhv data
        % Each profile goes into a row of LHV
        if (strcmp(Config.ProfPlots(iplot).Type, 'diff'))
          LHV(icase,:) = LHV_ALLZ(iz1:iz2,:) - LHV_CONTROL_ALLZ(iz1:iz2,:);
        else
          LHV(icase,:) = LHV_ALLZ(iz1:iz2,:);
        end
      end
    end

    % Latent heating rate is K/5min (RAMS files are saved every
    % 5 min of sim time), so need to multiply numbers the
    % resulting average profile by 12 to get to K/hr.
    LHV = LHV .* 12;
    Z = Z ./ 1000; % km

    fprintf('Writing plot file: %s\n', OutFile);
    PlotProfSet(X, Z, LHV, Xlabel, Zlabel, Ptitle, Pmarkers, Lstyles, Lgscales, LegText, LegLoc, Ptype, OutFile);
    fprintf('\n');
end
