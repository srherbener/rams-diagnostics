function [ ] = GenTsPlots(ConfigFile)
% GenTsPlots generate time series plots

% Read the config file to get the structure of how the data is laid out in
% the file system.
[ Config ] = ReadConfig(ConfigFile);
    
Pname = Config.Pexp.Ename;
Tstart = Config.Pexp.Tstart;
Tend = Config.Pexp.Tend;

Tdir = Config.TsavgDir;
Pdir = Config.PlotDir;

Flen = 5;
Times = (Tstart:Tend);
Tlen = length(Times);

Lcolors = { 'k' 'm' 'b' 'c' 'g' 'y' 'r' };

% Find and replace underscores in Ptitle, Ylabel with blank spaces
for iplot = 1:length(Config.TsavgPlots)
    Var = Config.TsavgPlots(iplot).Var;

    Ptitle = regexprep(sprintf('%s: %s', Pname, Config.TsavgPlots(iplot).Title), '_', ' ');
    LegLoc = Config.TsavgPlots(iplot).LegLoc;
    Ylim = [ Config.TsavgPlots(iplot).Ymin Config.TsavgPlots(iplot).Ymax ];
    Ylabel = regexprep(sprintf('%s (%s)', Config.TsavgPlots(iplot).Name, Config.TsavgPlots(iplot).Units), '_', ' ');
    OutFile = sprintf('%s/TS_%s.jpg', Pdir, Var);
    
    % make sure output directory exists
    if (exist(Pdir, 'dir') ~= 7)
        mkdir(Pdir);
    end
    
    for icase = 1:length(Config.Cases)
        Case = Config.Cases(icase).Cname;

        Hfile = sprintf('%s/%s_%s.h5', Tdir, Var, Case);
        fprintf('Reading HDF5 file: %s\n', Hfile);
        [ Hvar, Rcoords, Zcoords, Tcoords ] = ReadAzavgVar(Hfile, Var);
        TS = squeeze(Hvar);
        
        % smooth with a running mean of length 'Flen'
        [ TsAll(icase,:) ] = SmoothFillTseries(TS, Tlen, Flen);
        
        LegText(icase) = { Config.Cases(icase).Pname };
    end
    
    PlotTseriesSet( Times, TsAll, Ptitle, Ylim, Ylabel, Lcolors, LegText, LegLoc, OutFile );
    fprintf('\n');
end
