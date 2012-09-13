function [ ] = GenTsPlots(ConfigFile)
% GenTsPlots generate time series plots

% Read the config file to get the structure of how the data is laid out in
% the file system.
[ Config ] = ReadConfig(ConfigFile);
    
Pname   = Config.Pexp.Ename;
Ntsteps = Config.Pexp.Ntsteps;
Tstart  = Config.Pexp.Tstart;
Tinc    = Config.Pexp.Tinc;
Tunits  = Config.Pexp.Tunits;

Tdir = Config.TsavgDir;
Pdir = Config.PlotDir;

% For smoothing, length of a running mean
Flen = 5;

% Generate the time values
Times = zeros(1,Ntsteps);
for i = 1:Ntsteps
  Times(i) = Tstart + ((i-1)*Tinc);
end

Lcolors = { 'k' 'm' 'b' 'c' 'g' 'y' 'r' };

% Find and replace underscores in Ptitle, Ylabel with blank spaces
for iplot = 1:length(Config.TsavgPlots)
    Var = Config.TsavgPlots(iplot).Var;

    Ptitle = sprintf('%s: %s', Pname, Config.TsavgPlots(iplot).Title);
    LegLoc = Config.TsavgPlots(iplot).LegLoc;
    Ylim = [ Config.TsavgPlots(iplot).Ymin Config.TsavgPlots(iplot).Ymax ];
    Ylabel = sprintf('%s (%s)', Config.TsavgPlots(iplot).Name, Config.TsavgPlots(iplot).Units);
    OutFile = sprintf('%s/%s', Pdir, Config.TsavgPlots(iplot).OutFile);
    
    % make sure output directory exists
    if (exist(Pdir, 'dir') ~= 7)
        mkdir(Pdir);
    end

    ips = Config.TsavgPlots(iplot).PSnum;
    if (ips == 0)
      fprintf('WARNING: skipping TsavgPlot number %d due to no associated PlotSet\n', iplot)
    else
      for icase = 1:Config.PlotSets(ips).Ncases
        Case = Config.PlotSets(ips).Cases(icase).Cname;
        LegText(icase) = { Config.PlotSets(ips).Cases(icase).Legend };

        Hfile = sprintf('%s/%s_%s.h5', Tdir, Var, Case);
        fprintf('Reading HDF5 file: %s\n', Hfile);
        [ Hvar, Rcoords, Zcoords, Tcoords ] = ReadAzavgVar(Hfile, Var);
        TS = squeeze(Hvar);
        
        % smooth with a running mean of length 'Flen'
        [ TsAll(icase,:) ] = SmoothFillTseries(TS, Ntsteps, Flen);
      end
    end
    
    PlotTseriesSet( Times, TsAll, Ptitle, Ylim, Ylabel, Tunits, Lcolors, LegText, LegLoc, OutFile );
    fprintf('\n');
end
