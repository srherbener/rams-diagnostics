function [ ] = GenCloudStructPlots(ConfigFile)
% GenCloudStrucPlots function to plot cloud top temp vs cloud fraction

[ Config ] = ReadConfig(ConfigFile);

Pname   = Config.ExpName;

Tdir = Config.TsavgDir;
Pdir = Config.PlotDir;

% For smoothing, length of a running mean
Flen = 5;

Tstart = 120;
Tend = 433;
Tlen = (Tend - Tstart) + 1;

AxisProps(1).Name = 'FontSize';
AxisProps(1).Val = 20; 
AxisProps(2).Name = 'XDir';
AxisProps(2).Val = 'Reverse';
AxisProps(3).Name = 'Ylim';
AxisProps(3).Val = [ 0 1 ];
AxisProps(3).Name = 'Xlim';
AxisProps(3).Val = [ 0 15 ];

% make the plots
for iplot = 1:length(Config.TwoDimPlots)
    clear CtopAll;
    clear CfracAll;
    clear LegText;
    clear Lspecs;

    CtopVar = Config.TwoDimPlots(iplot).Xvar;
    CfracVar = Config.TwoDimPlots(iplot).Yvar;

    Ptitle = sprintf('%s: %s', Pname, Config.TwoDimPlots(iplot).Title);
    Xlabel = Config.TwoDimPlots(iplot).Xlabel;
    Ylabel = Config.TwoDimPlots(iplot).Ylabel;
    LegLoc = Config.TwoDimPlots(iplot).LegLoc;
    OutFile = sprintf('%s/%s', Pdir, Config.TwoDimPlots(iplot).OutFile);
    
    % make sure output directory exists
    if (exist(Pdir, 'dir') ~= 7)
        mkdir(Pdir);
    end

    ips = Config.TwoDimPlots(iplot).PSnum;
    if (ips == 0)
      fprintf('WARNING: skipping TwoDimPlot number %d due to no associated PlotSet\n', iplot)
    else
      for icase = 1:Config.PlotSets(ips).Ncases
        Case = Config.PlotSets(ips).Cases(icase).Cname;
        LegText(icase) = { Config.PlotSets(ips).Cases(icase).Legend };
        Lspecs(icase) = { Config.PlotSets(ips).Cases(icase).Lspec };

        CtopFile = sprintf('%s/%s_%s.h5', Tdir, CtopVar, Case);
        fprintf('Reading HDF5 file: %s\n', CtopFile);
        CTOP = squeeze(hdf5read(CtopFile, CtopVar));
        CTOP = CTOP(Tstart:Tend);
        [ CtopAll(icase,:) ] = SmoothFillTseries(CTOP, Tlen, Flen);
    
        CfracFile = sprintf('%s/%s_%s.h5', Tdir, CfracVar, Case);
        fprintf('Reading HDF5 file: %s\n', CfracFile);
        CFRAC = squeeze(hdf5read(CfracFile, CfracVar));
        CFRAC = CFRAC(Tstart:Tend);
        [ CfracAll(icase,:) ] = SmoothFillTseries(CFRAC, Tlen, Flen);
      end
    end
    
    fprintf('Writing file: %s\n', OutFile);
    Plot2dSet( CtopAll, CfracAll, Ptitle, Xlabel, Ylabel, Lspecs, LegText, LegLoc, AxisProps, OutFile );
    fprintf('\n');
end

end
