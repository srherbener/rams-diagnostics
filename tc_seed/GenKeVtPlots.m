function [ ] = GenKeVtPlots(ConfigFile)
% GenKeVtPlots function to plot total KE vs max azimuthally averaged tangential wind

[ Config ] = ReadConfig(ConfigFile);

UndefVal = Config.UndefVal;

Pname   = Config.ExpName;

Tdir = Config.TsavgDir;
Pdir = Config.PlotDir;

% For smoothing, length of a running mean
Flen = 5;

Tstart = 1;
Tend = Config.TsPlotSpecs.Ntsteps;
Tlen = (Tend - Tstart) + 1;

CntlTstart = Config.TsPlotSpecs.ControlStart;
CntlTend = CntlTstart + (Config.TsPlotSpecs.Ntsteps - 1);

% For plotting
Lcolors = { 'k' 'm' 'b' 'c' 'g' 'y' 'r' };

AxisProps(1).Name = 'FontSize';
AxisProps(1).Val = 20; 

% make the plots
for iplot = 1:length(Config.TwoDimPlots)
    clear VtAll;
    clear KeAll;
    clear LegText;

    VtVar = Config.TwoDimPlots(iplot).Xvar;
    KeVar = Config.TwoDimPlots(iplot).Yvar;

    Ptitle = sprintf('%s: %s', Pname, Config.TwoDimPlots(iplot).Title);
    Xlabel = Config.TwoDimPlots(iplot).Xlabel;
    Xscale = Config.TwoDimPlots(iplot).Xscale;
    Ylabel = Config.TwoDimPlots(iplot).Ylabel;
    Yscale = Config.TwoDimPlots(iplot).Yscale;
    LegLoc = Config.TwoDimPlots(iplot).LegLoc;
    OutFile = sprintf('%s/%s', Pdir, Config.TwoDimPlots(iplot).OutFile);
    
    % make sure output directory exists
    if (exist(Pdir, 'dir') ~= 7)
        mkdir(Pdir);
    end

    ihist = 0;
    ips = Config.TwoDimPlots(iplot).PSnum;
    if (ips == 0)
      fprintf('WARNING: skipping TwoDimPlot number %d due to no associated PlotSet\n', iplot)
    else
      for icase = 1:Config.PlotSets(ips).Ncases
        Case = Config.PlotSets(ips).Cases(icase).Cname;
        LegText(icase) = { Config.PlotSets(ips).Cases(icase).Legend };

        VtFile = sprintf('%s/%s_%s.h5', Tdir, VtVar, Case);
        fprintf('Reading HDF5 file: %s\n', VtFile);
        VT = squeeze(hdf5read(VtFile, VtVar));
        if (strcmp(Case, Config.ControlCase))
          VT = VT(CntlTstart:CntlTend);
        else
          VT = VT(Tstart:Tend);
        end
        VT = VT * Xscale;
        [ VtAll(icase,:) ] = SmoothFillTseries(VT, Tlen, Flen);
    
        KeFile = sprintf('%s/%s_%s.h5', Tdir, KeVar, Case);
        fprintf('Reading HDF5 file: %s\n', KeFile);
        KE = squeeze(hdf5read(KeFile, KeVar));
        if (strcmp(Case, Config.ControlCase))
          KE = KE(CntlTstart:CntlTend);
        else
          KE = KE(Tstart:Tend);
        end
        KE = KE * Yscale;
        [ KeAll(icase,:) ] = SmoothFillTseries(KE, Tlen, Flen);
      end
    end

    Plot2dSet( VtAll, KeAll, Ptitle, Xlabel, Ylabel, Lcolors, LegText, LegLoc, AxisProps, OutFile );
    fprintf('\n');
end

end
