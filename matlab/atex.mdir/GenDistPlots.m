function [ ] = GenDistPlots(ConfigFile)
% GenDistPlots generate histogram distribution plots

% Read the config file to get the structure of how the data is laid out in
% the file system.
[ Config ] = ReadConfig(ConfigFile);
    
Pname   = Config.ExpName;

Ddir = Config.DiagDir;
Pdir = Config.PlotDir;

Lcolors = { 'k' 'm' 'b' 'c' 'g' 'y' 'r' };

% Find and replace underscores in Ptitle, Ylabel with blank spaces
for iplot = 1:length(Config.DistPlots)
    clear H;
    clear B;
    clear LegText;

    Var = Config.DistPlots(iplot).Var;

    Ptitle = sprintf('%s: %s', Pname, Config.DistPlots(iplot).Title);
    Xlabel = Config.DistPlots(iplot).Xlabel;
    Ylabel = Config.DistPlots(iplot).Ylabel;
    LegLoc = Config.DistPlots(iplot).LegLoc;
    OutFile = sprintf('%s/%s', Pdir, Config.DistPlots(iplot).OutFile);
    
    % make sure output directory exists
    if (exist(Pdir, 'dir') ~= 7)
        mkdir(Pdir);
    end

    ihist = 0;
    ips = Config.DistPlots(iplot).PSnum;
    if (ips == 0)
      fprintf('WARNING: skipping DistPlot number %d due to no associated PlotSet\n', iplot)
    else
      for icase = 1:Config.PlotSets(ips).Ncases
        Case = Config.PlotSets(ips).Cases(icase).Cname;
        LegText(icase) = { Config.PlotSets(ips).Cases(icase).Legend };
        Ylim = [ Config.DistPlots(iplot).Ymin Config.DistPlots(iplot).Ymax ];

        Hfile = sprintf('%s/%s_%s.h5', Ddir, Var, Case);
        fprintf('Reading HDF5 file: %s\n', Hfile);

        % histograms are in the columns
        ihist = ihist + 1;
        H(:,ihist) = hdf5read(Hfile, 'AvgDist');
        B(:,ihist) = hdf5read(Hfile, 'Bins');
      end
    end

    % just look at the lower rain rates --> 0 to 2.9 mm/hr
    Hists = H(1:29,:);
    Bins  = B(1:29,:);

    PlotDistSet( Hists, Bins, Ptitle, Xlabel, Ylabel, Ylim, Lcolors, LegText, LegLoc, OutFile );
    fprintf('\n');
end
