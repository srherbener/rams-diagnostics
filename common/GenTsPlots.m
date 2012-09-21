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

% BT is an array with numeric entries for:
%    1 - year
%    2 - month
%    3 - day
%    4 - hour
%    5 - minute
%    6 - second
BT = Config.Pexp.BaseTime;
BaseTime = datenum(BT(1), BT(2), BT(3), BT(4), BT(5), BT(6));
StartTime = datestr(BaseTime, 'mm/dd/yyyy HH:MM')

TsStart = Config.Pexp.TsStart;
TsPeriod = Config.Pexp.TsPeriod;

Tdir = Config.TsavgDir;
Pdir = Config.PlotDir;

% For smoothing, length of a running mean
Flen = 5;

% Generate the time values
Times = zeros(1,Ntsteps);
it = 0;
for i = 1:Ntsteps
  % Times will be in hours
  Times(i) = Tstart + ((i-1)*Tinc);

  if (mod(i-TsStart,TsPeriod) == 0)
    it = it + 1;
    % returns date and time strings
    [ Ds, Ts ] = TimeToString(BT(1), BT(2), BT(3), BT(4), BT(5), BT(6), Times(i));

    Tticks(it) = Times(i);
    Tlabels{it} =  Ts;
  end
end

Lcolors = { 'k' 'm' 'b' 'c' 'g' 'y' 'r' };

% Find and replace underscores in Ptitle, Ylabel with blank spaces
for iplot = 1:length(Config.TsavgPlots)
    clear TsAll;
    clear LegText;

    Var = Config.TsavgPlots(iplot).Var;

    % If doing a diff plot, read in the control profile
    if (strcmp(Config.TsavgPlots(iplot).Type, 'diff'))
      Case = Config.ControlCase;
      Hfile = sprintf('%s/%s_%s.h5', Tdir, Var, Case);
      fprintf('Reading Control Case: %s\n', Case);
      fprintf('  HDF5 file: %s\n', Hfile);
      [ Hvar, Rcoords, Zcoords, Tcoords ] = ReadAzavgVar(Hfile, Var);
      TS_CNTL = squeeze(Hvar);
    end 

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
        
        % If doing a diff type plot, subtract off the control values
        if (strcmp(Config.TsavgPlots(iplot).Type, 'diff'))
          TS = TS - TS_CNTL;
        end

        % smooth with a running mean of length 'Flen'
        [ TsAll(icase,:) ] = SmoothFillTseries(TS, Ntsteps, Flen);
      end
    end
    
    PlotTseriesSet( Times, TsAll, Ptitle, Ylim, Ylabel, StartTime, Tticks, Tlabels, Tunits, Lcolors, LegText, LegLoc, OutFile );
    fprintf('\n');
end
