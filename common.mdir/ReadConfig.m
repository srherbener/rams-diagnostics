function [ Cdata ] = ReadConfig ( Cfile )
% ReadConfig read diagnostic configuration file for a set of simulations
%
% This function will configuration file for processing
% a set of experimental simulation runs.
%
% The config file has the format:
%    Case: <directory>
%    TimeDir: <directory>
%
% Case: and TimeDir: define all of the directories where hdf5 files live
% 
% Case: defines separate experiments (simulations)
%
% TimeDir: defines multiple directories under a single experiment where
%       hdf5 files live
%
% If you have
%   Case: RUN1
%   Case: RUN2
%   Case: RUN3
%   TimeDir: T1
%   TimeDir: T2
%
% then the directories where the hdf5 files live are:
%
%   For RUN1:
%     RUN1/HDF5/T1
%     RUN1/HDF5/T2
%
%   For RUN2:
%     RUN2/HDF5/T1
%     RUN2/HDF5/T2
%
%   For RUN3:
%     RUN3/HDF5/T1
%     RUN3/HDF5/T2
%

% Read in entire file and place into a character array.
fid = fopen(char(Cfile));
InFile = textscan(fid, '%s', 'Delimiter', '\n'); % read in the entire line
fclose(fid);

% Convert to char array, entire file is in first element of InFile.
InLines = char(InFile{1});

% Grab keywords and their data and load into structs
i_case = 0;
i_tdir = 0;
i_aeof = 0;
i_hmeas = 0;
i_pmeas = 0;
i_smeas = 0;
i_smeas1d = 0;
i_hist2d = 0;
i_counts = 0;
i_tavg = 0;
i_dpts = 0;
i_vint_ts = 0;
i_aeofplot = 0;
i_tsplot = 0;
i_dplot = 0;
i_lplot = 0;
i_hmp3d = 0;
i_hmslice = 0;
i_pplot = 0;
i_splot = 0;
i_splot1d = 0;
i_pcplot = 0;
i_pts_plot = 0;
i_prs_plot = 0;
i_pset = 0;
i_pvar = 0;
i_pds = 0;
for i = 1:size(InLines,1)
  % Convert a line to a list of space separated fields
  [ Fields ] = Line2Fields(InLines(i,:), ' ');

  switch Fields{1}
    case 'AzavgDir:'
      Cdata.AzavgDir = Fields{2};
    case 'TsavgDir:'
      Cdata.TsavgDir = Fields{2};
    case 'DiagDir:'
      Cdata.DiagDir = Fields{2};
    case 'PlotDir:'
      Cdata.PlotDir = Fields{2};
    case 'FigDir:'
      Cdata.FigDir = Fields{2};
    case 'EofDir:'
      Cdata.EofDir = Fields{2};
    case 'ControlCase:'
      Cdata.ControlCase = Fields{2};
    case 'SpinUpCase:'
      Cdata.SpinUpCase = Fields{2};
    case 'ExpName:'
      Cdata.ExpName = Fields{2};
    case 'UndefVal:'
      Cdata.UndefVal = sscanf(Fields{2}, '%f');
    case 'SstVal:'
      Cdata.SstVal = sscanf(Fields{2}, '%f');
    case 'Case:'
      i_case = i_case + 1;
      Cdata.Cases(i_case).Cname = Fields{2};
      Cdata.Cases(i_case).Pname = Fields{3};
    case 'TimeDir:'
      i_tdir = i_tdir + 1;
      Cdata.Tdirs{i_tdir} = Fields{2};
    case 'AzavgEofConfig:'
      Cdata.AzavgEofConfig.NumEv    = sscanf(Fields{2}, '%d');
      Cdata.AzavgEofConfig.Nstar    = sscanf(Fields{3}, '%d');
    case 'AzavgEof:'
      i_aeof = i_aeof + 1;
      Cdata.AzavgEof(i_aeof).Name    = Fields{2};
      Cdata.AzavgEof(i_aeof).Fprefix = Fields{3};
      Cdata.AzavgEof(i_aeof).Rvar    = Fields{4};
      Cdata.AzavgEof(i_aeof).Rmin    = sscanf(Fields{5}, '%f');
      Cdata.AzavgEof(i_aeof).Rmax    = sscanf(Fields{6}, '%f');
      Cdata.AzavgEof(i_aeof).Zmin    = sscanf(Fields{7}, '%f');
      Cdata.AzavgEof(i_aeof).Zmax    = sscanf(Fields{8}, '%f');
      Cdata.AzavgEof(i_aeof).Tmin    = sscanf(Fields{9}, '%f');
      Cdata.AzavgEof(i_aeof).Tmax    = sscanf(Fields{10}, '%f');
    case 'Hmeas:'
      i_hmeas = i_hmeas + 1;
      Cdata.Hmeas(i_hmeas).Name    = Fields{2};
      Cdata.Hmeas(i_hmeas).InDir   = Fields{3};
      Cdata.Hmeas(i_hmeas).Fprefix = Fields{4};
      Cdata.Hmeas(i_hmeas).Rvar    = Fields{5};
      Cdata.Hmeas(i_hmeas).Method  = Fields{6};
    case 'Pmeas:'
      i_pmeas = i_pmeas + 1;
      Cdata.Pmeas(i_pmeas).Name    = Fields{2};
      Cdata.Pmeas(i_pmeas).InDir   = Fields{3};
      Cdata.Pmeas(i_pmeas).Fprefix = Fields{4};
      Cdata.Pmeas(i_pmeas).Rvar    = Fields{5};
      Cdata.Pmeas(i_pmeas).Method  = Fields{6};
      Cdata.Pmeas(i_pmeas).Rmin    = sscanf(Fields{7},  '%f');
      Cdata.Pmeas(i_pmeas).Rmax    = sscanf(Fields{8},  '%f');
      Cdata.Pmeas(i_pmeas).Tmin    = sscanf(Fields{9},  '%f');
      Cdata.Pmeas(i_pmeas).Tmax    = sscanf(Fields{10}, '%f');
    case 'Smeas:'
      i_smeas = i_smeas + 1;
      Cdata.Smeas(i_smeas).PSname  = Fields{2};
      Cdata.Smeas(i_smeas).Name    = Fields{3};
      Cdata.Smeas(i_smeas).InDir   = Fields{4};
      Cdata.Smeas(i_smeas).Fprefix = Fields{5};
      Cdata.Smeas(i_smeas).Rvar    = Fields{6};
      Cdata.Smeas(i_smeas).Xmin    = sscanf(Fields{7}, '%f');
      Cdata.Smeas(i_smeas).Xmax    = sscanf(Fields{8}, '%f');
      Cdata.Smeas(i_smeas).Xgroup  = sscanf(Fields{9}, '%d');
      Cdata.Smeas(i_smeas).Ymin    = sscanf(Fields{10}, '%f');
      Cdata.Smeas(i_smeas).Ymax    = sscanf(Fields{11}, '%f');
      Cdata.Smeas(i_smeas).Ygroup  = sscanf(Fields{12}, '%d');
      ix = 1;
      for j = 13:length(Fields)  % grab the rest of the fields
        Cdata.Smeas(i_smeas).Xvals(ix) = sscanf(Fields{j}, '%f');
        ix = ix + 1;
      end

      Cdata.Smeas(i_smeas).PSnum   =  -1;
    case 'Smeas1d:'
      i_smeas1d = i_smeas1d + 1;
      Cdata.Smeas1d(i_smeas1d).PSname  = Fields{2};
      Cdata.Smeas1d(i_smeas1d).Name    = Fields{3};
      Cdata.Smeas1d(i_smeas1d).InDir   = Fields{4};
      Cdata.Smeas1d(i_smeas1d).Fprefix = Fields{5};
      Cdata.Smeas1d(i_smeas1d).Rvar    = Fields{6};
      Cdata.Smeas1d(i_smeas1d).Xmin    = sscanf(Fields{7}, '%f');
      Cdata.Smeas1d(i_smeas1d).Xmax    = sscanf(Fields{8}, '%f');
      Cdata.Smeas1d(i_smeas1d).Xgroup  = sscanf(Fields{9}, '%d');
      ix = 1;
      for j = 10:length(Fields)  % grab the rest of the fields
        Cdata.Smeas1d(i_smeas1d).Xvals(ix) = sscanf(Fields{j}, '%f');
        ix = ix + 1;
      end

      Cdata.Smeas1d(i_smeas1d).PSnum   =  -1;
    case 'Hist2d:'
      i_hist2d = i_hist2d + 1;
      Cdata.Hist2d(i_hist2d).Name    = Fields{2};
      Cdata.Hist2d(i_hist2d).Fprefix = Fields{3};
      Cdata.Hist2d(i_hist2d).Var     = Fields{4};
      Cdata.Hist2d(i_hist2d).Xmin    = sscanf(Fields{5}, '%f');
      Cdata.Hist2d(i_hist2d).Xmax    = sscanf(Fields{6}, '%f');
      Cdata.Hist2d(i_hist2d).Xgroup  = sscanf(Fields{7}, '%d');
      Cdata.Hist2d(i_hist2d).Ymin    = sscanf(Fields{8}, '%f');
      Cdata.Hist2d(i_hist2d).Ymax    = sscanf(Fields{9}, '%f');
      Cdata.Hist2d(i_hist2d).Ygroup  = sscanf(Fields{10}, '%d');
      Cdata.Hist2d(i_hist2d).Tmin    = sscanf(Fields{11}, '%f');
      Cdata.Hist2d(i_hist2d).Tmax    = sscanf(Fields{12}, '%f');
    case 'CompCounts:'
      i_counts = i_counts + 1;
      Cdata.CompCounts(i_counts).PSname  = Fields{2};
      Cdata.CompCounts(i_counts).Name    = Fields{3};
      Cdata.CompCounts(i_counts).InDir   = Fields{4};
      Cdata.CompCounts(i_counts).Fprefix = Fields{5};
      Cdata.CompCounts(i_counts).Rvar    = Fields{6};

      Cdata.CompCounts(i_counts).PSnum   =  -1;
    case 'Tavg:'
      i_tavg = i_tavg + 1;
      Cdata.Tavg(i_tavg).Name    = Fields{2};
      Cdata.Tavg(i_tavg).InDir   = Fields{3};
      Cdata.Tavg(i_tavg).Fprefix = Fields{4};
      Cdata.Tavg(i_tavg).Rvar    = Fields{5};
      Cdata.Tavg(i_tavg).Tmin    = sscanf(Fields{6}, '%f');
      Cdata.Tavg(i_tavg).Tmax    = sscanf(Fields{7}, '%f');
    case 'DeltapTs:'
      i_dpts = i_dpts + 1;
      Cdata.DeltapTs(i_dpts).Name    = Fields{2};
      Cdata.DeltapTs(i_dpts).InDir   = Fields{3};
      Cdata.DeltapTs(i_dpts).Fprefix = Fields{4};
      Cdata.DeltapTs(i_dpts).Rvar    = Fields{5};
      Cdata.DeltapTs(i_dpts).Rmin    = sscanf(Fields{6}, '%f');
      Cdata.DeltapTs(i_dpts).Rmax    = sscanf(Fields{7}, '%f');
    case 'VintTs:'
      i_vint_ts = i_vint_ts + 1;
      Cdata.VintTs(i_vint_ts).Name    = Fields{2};
      Cdata.VintTs(i_vint_ts).InDir   = Fields{3};
      Cdata.VintTs(i_vint_ts).Fprefix = Fields{4};
      Cdata.VintTs(i_vint_ts).Rvar    = Fields{5};
    case 'TsPlotSpecs:'
      Cdata.TsPlotSpecs.Ntsteps      = sscanf(Fields{2}, '%d');
      Cdata.TsPlotSpecs.Tstart       = sscanf(Fields{3}, '%f');
      Cdata.TsPlotSpecs.Tinc         = sscanf(Fields{4}, '%f');
      Cdata.TsPlotSpecs.Tunits       = Fields{5};
      Cdata.TsPlotSpecs.ControlStart = sscanf(Fields{6}, '%d');
      Cdata.TsPlotSpecs.BaseTime     = sscanf(Fields{7}, '%d:%d:%d:%d:%d:%d', 6);
      Cdata.TsPlotSpecs.TsStart      = sscanf(Fields{8}, '%d');
      Cdata.TsPlotSpecs.TsPeriod     = sscanf(Fields{9}, '%d');
    case 'PlotSet:'
      i_pset = i_pset + 1;
      Cdata.PlotSets(i_pset).Name   = Fields{2};
      Cdata.PlotSets(i_pset).Ncases = sscanf(Fields{3}, '%d');
      j = 4; % next field
      for ips = 1:Cdata.PlotSets(i_pset).Ncases
        Cdata.PlotSets(i_pset).Cases(ips).Cname   = Fields{j};
        Cdata.PlotSets(i_pset).Cases(ips).Legend  = regexprep(Fields{j+1}, '_', ' ');
        Cdata.PlotSets(i_pset).Cases(ips).Lstyle  = Fields{j+2};
        Cdata.PlotSets(i_pset).Cases(ips).Lgscale = sscanf(Fields{j+3}, '%f');
        j = j + 4;
      end
    case 'PlotVar:'
      i_pvar = i_pvar + 1;
      Cdata.PlotVars(i_pvar).Name    = Fields{2};
      Cdata.PlotVars(i_pvar).Var     = Fields{3};
      Cdata.PlotVars(i_pvar).Fprefix = Fields{4};
      Cdata.PlotVars(i_pvar).Label   = regexprep(Fields{5}, '_', ' ');
      Cdata.PlotVars(i_pvar).Units   = regexprep(Fields{6}, '_', ' ');
      Cdata.PlotVars(i_pvar).Min     = sscanf(Fields{7}, '%f');
      Cdata.PlotVars(i_pvar).Max     = sscanf(Fields{8}, '%f');
      Cdata.PlotVars(i_pvar).Scale   = sscanf(Fields{9}, '%f');
    case 'PlotDselect:'
      i_pds = i_pds + 1;
      Cdata.PlotDselects(i_pds).Name = Fields{2};
      Cdata.PlotDselects(i_pds).Xmin = sscanf(Fields{3}, '%f');
      Cdata.PlotDselects(i_pds).Xmax = sscanf(Fields{4}, '%f');
      Cdata.PlotDselects(i_pds).Ymin = sscanf(Fields{5}, '%f');
      Cdata.PlotDselects(i_pds).Ymax = sscanf(Fields{6}, '%f');
      Cdata.PlotDselects(i_pds).Zmin = sscanf(Fields{7}, '%f');
      Cdata.PlotDselects(i_pds).Zmax = sscanf(Fields{8}, '%f');
      Cdata.PlotDselects(i_pds).Tmin = sscanf(Fields{9}, '%f');
      Cdata.PlotDselects(i_pds).Tmax = sscanf(Fields{10}, '%f');
    case 'TsavgPlot:'
      i_tsplot = i_tsplot + 1;
      Cdata.TsavgPlots(i_tsplot).PSname  = Fields{2};
      Cdata.TsavgPlots(i_tsplot).Fprefix = Fields{3};
      Cdata.TsavgPlots(i_tsplot).Var     = Fields{4};
      Cdata.TsavgPlots(i_tsplot).Ptype   = Fields{5};
      Cdata.TsavgPlots(i_tsplot).Ttype   = Fields{6};
      Cdata.TsavgPlots(i_tsplot).Name    = regexprep(Fields{7}, '_', ' ');
      Cdata.TsavgPlots(i_tsplot).Units   = regexprep(Fields{8}, '_', ' ');
      Cdata.TsavgPlots(i_tsplot).Title   = ParseTitle(Fields{9});
      Cdata.TsavgPlots(i_tsplot).LegLoc  = Fields{10};
      Cdata.TsavgPlots(i_tsplot).Ymin    = sscanf(Fields{11},  '%f');
      Cdata.TsavgPlots(i_tsplot).Ymax    = sscanf(Fields{12},  '%f');
      Cdata.TsavgPlots(i_tsplot).OutFile = Fields{13};
      % Need the following to get the assignment from AssociateStructs to work
      Cdata.TsavgPlots(i_tsplot).PSnum   = -1;
    case 'DistPlot:'
      i_dplot = i_dplot + 1;
      Cdata.DistPlots(i_dplot).PSname  = Fields{2};
      Cdata.DistPlots(i_dplot).Var     = Fields{3};
      Cdata.DistPlots(i_dplot).Title   = ParseTitle(Fields{4});
      Cdata.DistPlots(i_dplot).Xlabel  = regexprep(Fields{5}, '_', ' ');
      Cdata.DistPlots(i_dplot).Ylabel  = regexprep(Fields{6}, '_', ' ');
      Cdata.DistPlots(i_dplot).LegLoc  = Fields{7};
      Cdata.DistPlots(i_dplot).Ymin    = sscanf(Fields{8},  '%f');
      Cdata.DistPlots(i_dplot).Ymax    = sscanf(Fields{9},  '%f');
      Cdata.DistPlots(i_dplot).OutFile = Fields{10};

      Cdata.DistPlots(i_dplot).PSnum   = -1;
    case 'LinePlot:'
      i_lplot = i_lplot + 1;
      Cdata.LinePlots(i_lplot).PSname  = Fields{2};
      Cdata.LinePlots(i_lplot).XVname  = Fields{3};
      Cdata.LinePlots(i_lplot).YVname  = Fields{4};
      Cdata.LinePlots(i_lplot).DSname  = Fields{5};
      Cdata.LinePlots(i_lplot).Smooth  = Fields{6};
      Cdata.LinePlots(i_lplot).Title   = ParseTitle(Fields{7});
      Cdata.LinePlots(i_lplot).LegLoc  = Fields{8};
      Cdata.LinePlots(i_lplot).AddMeas = Fields{9};
      Cdata.LinePlots(i_lplot).OutFile = Fields{10};

      Cdata.LinePlots(i_lplot).PSnum   = -1;
      Cdata.LinePlots(i_lplot).XVnum   = -1;
      Cdata.LinePlots(i_lplot).YVnum   = -1;
      Cdata.LinePlots(i_lplot).DSnum   = -1;
    case 'HmeasPlot3d:'
      i_hmp3d = i_hmp3d + 1;
      Cdata.HmeasPlot3d(i_hmp3d).Name    = Fields{2};
      Cdata.HmeasPlot3d(i_hmp3d).Var     = Fields{3};
      Cdata.HmeasPlot3d(i_hmp3d).Fprefix = Fields{4};
      Cdata.HmeasPlot3d(i_hmp3d).Descrip = regexprep(Fields{5}, '_', ' ');
      Cdata.HmeasPlot3d(i_hmp3d).Units   = Fields{6};
      Cdata.HmeasPlot3d(i_hmp3d).Ptype   = Fields{7};
      Cdata.HmeasPlot3d(i_hmp3d).Vtype   = Fields{8};
      Cdata.HmeasPlot3d(i_hmp3d).Isurf   = sscanf(Fields{9}, '%f');
      Cdata.HmeasPlot3d(i_hmp3d).Rmin    = sscanf(Fields{10},'%f');
      Cdata.HmeasPlot3d(i_hmp3d).Rmax    = sscanf(Fields{11},'%f');
      Cdata.HmeasPlot3d(i_hmp3d).Zmin    = sscanf(Fields{12},'%f');
      Cdata.HmeasPlot3d(i_hmp3d).Zmax    = sscanf(Fields{13},'%f');
      Cdata.HmeasPlot3d(i_hmp3d).Tmin    = sscanf(Fields{14},'%f');
      Cdata.HmeasPlot3d(i_hmp3d).Tmax    = sscanf(Fields{15},'%f');
      Cdata.HmeasPlot3d(i_hmp3d).Cmin    = sscanf(Fields{16},'%f');
      Cdata.HmeasPlot3d(i_hmp3d).Cmax    = sscanf(Fields{17},'%f');
      Cdata.HmeasPlot3d(i_hmp3d).Flevel  = sscanf(Fields{18},'%d');
    case 'HmeasSlicePlot:'
      i_hmslice = i_hmslice + 1;
      Cdata.HmeasSlicePlots(i_hmslice).Name    = Fields{2};
      Cdata.HmeasSlicePlots(i_hmslice).PSname  = Fields{3};
      Cdata.HmeasSlicePlots(i_hmslice).Var     = Fields{4};
      Cdata.HmeasSlicePlots(i_hmslice).Fprefix = Fields{5};
      Cdata.HmeasSlicePlots(i_hmslice).Descrip = regexprep(Fields{6}, '_', ' ');
      Cdata.HmeasSlicePlots(i_hmslice).Units   = Fields{7};
      Cdata.HmeasSlicePlots(i_hmslice).Ptype   = Fields{8};
      Cdata.HmeasSlicePlots(i_hmslice).Vtype   = Fields{9};
      Cdata.HmeasSlicePlots(i_hmslice).Stype   = Fields{10};
      Cdata.HmeasSlicePlots(i_hmslice).S1      = sscanf(Fields{11},'%f');
      Cdata.HmeasSlicePlots(i_hmslice).S2      = sscanf(Fields{12},'%f');
      Cdata.HmeasSlicePlots(i_hmslice).S3      = sscanf(Fields{13},'%f');
      Cdata.HmeasSlicePlots(i_hmslice).Rmin    = sscanf(Fields{14},'%f');
      Cdata.HmeasSlicePlots(i_hmslice).Rmax    = sscanf(Fields{15},'%f');
      Cdata.HmeasSlicePlots(i_hmslice).Zmin    = sscanf(Fields{16},'%f');
      Cdata.HmeasSlicePlots(i_hmslice).Zmax    = sscanf(Fields{17},'%f');
      Cdata.HmeasSlicePlots(i_hmslice).Tmin    = sscanf(Fields{18},'%f');
      Cdata.HmeasSlicePlots(i_hmslice).Tmax    = sscanf(Fields{19},'%f');
      Cdata.HmeasSlicePlots(i_hmslice).Cmin    = sscanf(Fields{20},'%f');
      Cdata.HmeasSlicePlots(i_hmslice).Cmax    = sscanf(Fields{21},'%f');

      Cdata.HmeasSlicePlots(i_hmslice).PSnum   = -1;
    case 'ProfPlot:'
      i_pplot = i_pplot + 1;
      Cdata.ProfPlots(i_pplot).PSname  = Fields{2};
      Cdata.ProfPlots(i_pplot).Fprefix = Fields{3};
      Cdata.ProfPlots(i_pplot).Var     = Fields{4};
      Cdata.ProfPlots(i_pplot).Type    = Fields{5};
      Cdata.ProfPlots(i_pplot).Title   = ParseTitle(Fields{6});
      Cdata.ProfPlots(i_pplot).Xlabel  = regexprep(Fields{7}, '_', ' ');
      Cdata.ProfPlots(i_pplot).Zlabel  = regexprep(Fields{8}, '_', ' ');
      Cdata.ProfPlots(i_pplot).LegLoc  = Fields{9};
      Cdata.ProfPlots(i_pplot).Xspec   = sscanf(Fields{10}, '%f:%f:%f', 3);
      Cdata.ProfPlots(i_pplot).Zmin    = sscanf(Fields{11}, '%f');
      Cdata.ProfPlots(i_pplot).Zmax    = sscanf(Fields{12}, '%f');
      Cdata.ProfPlots(i_pplot).OutFile = Fields{13};

      Cdata.ProfPlots(i_pplot).PSnum   = -1;
    case 'SlopePlot:'
      i_splot = i_splot + 1;
      Cdata.SlopePlots(i_splot).InFile  = Fields{2};
      Cdata.SlopePlots(i_splot).Svar    = Fields{3};
      Cdata.SlopePlots(i_splot).Cvar    = Fields{4};
      Cdata.SlopePlots(i_splot).XLvar   = Fields{5};
      Cdata.SlopePlots(i_splot).XUvar   = Fields{6};
      Cdata.SlopePlots(i_splot).YLvar   = Fields{7};
      Cdata.SlopePlots(i_splot).YUvar   = Fields{8};
      Cdata.SlopePlots(i_splot).Tvar    = Fields{9};
      Cdata.SlopePlots(i_splot).Nsamp   = sscanf(Fields{10}, '%d');
      Cdata.SlopePlots(i_splot).Tcor    = sscanf(Fields{11}, '%f');
      Cdata.SlopePlots(i_splot).Title   = ParseTitle(Fields{12});
      Cdata.SlopePlots(i_splot).Xlabel  = regexprep(Fields{13}, '_', ' ');
      Cdata.SlopePlots(i_splot).Ylabel  = regexprep(Fields{14}, '_', ' ');
      Cdata.SlopePlots(i_splot).Tmin    = sscanf(Fields{15}, '%f');
      Cdata.SlopePlots(i_splot).Tmax    = sscanf(Fields{16}, '%f');
      Cdata.SlopePlots(i_splot).Cmin    = sscanf(Fields{17}, '%f');
      Cdata.SlopePlots(i_splot).Cmax    = sscanf(Fields{18}, '%f');
      Cdata.SlopePlots(i_splot).OutFile = Fields{19};
    case 'Slope1dPlot:'
      i_splot1d = i_splot1d + 1;
      Cdata.Slope1dPlots(i_splot1d).PSname  = Fields{2};
      Cdata.Slope1dPlots(i_splot1d).InDir   = Fields{3};
      Cdata.Slope1dPlots(i_splot1d).Svar    = Fields{4};
      Cdata.Slope1dPlots(i_splot1d).Cvar    = Fields{5};
      Cdata.Slope1dPlots(i_splot1d).XLvar   = Fields{6};
      Cdata.Slope1dPlots(i_splot1d).XUvar   = Fields{7};
      Cdata.Slope1dPlots(i_splot1d).Tvar    = Fields{8};
      Cdata.Slope1dPlots(i_splot1d).Nsamp   = sscanf(Fields{9}, '%d');
      Cdata.Slope1dPlots(i_splot1d).Tcor    = sscanf(Fields{10}, '%f');
      Cdata.Slope1dPlots(i_splot1d).Title   = ParseTitle(Fields{11});
      Cdata.Slope1dPlots(i_splot1d).Xlabel  = regexprep(Fields{12}, '_', ' ');
      Cdata.Slope1dPlots(i_splot1d).Ylabel  = regexprep(Fields{13}, '_', ' ');
      Cdata.Slope1dPlots(i_splot1d).Xmin    = sscanf(Fields{14}, '%f');
      Cdata.Slope1dPlots(i_splot1d).Xmax    = sscanf(Fields{15}, '%f');
      Cdata.Slope1dPlots(i_splot1d).Tmin    = sscanf(Fields{16}, '%f');
      Cdata.Slope1dPlots(i_splot1d).Tmax    = sscanf(Fields{17}, '%f');
      Cdata.Slope1dPlots(i_splot1d).Cmin    = sscanf(Fields{18}, '%f');
      Cdata.Slope1dPlots(i_splot1d).Cmax    = sscanf(Fields{19}, '%f');
      Cdata.Slope1dPlots(i_splot1d).OutFile = Fields{20};
      
      Cdata.Slope1dPlots(i_splot1d).PSnum = -1;
    case 'PcolorPlot:'
      i_pcplot = i_pcplot + 1;
      Cdata.PcolorPlots(i_pcplot).Fprefix = Fields{2};
      Cdata.PcolorPlots(i_pcplot).Cmin    = sscanf(Fields{3}, '%f');
      Cdata.PcolorPlots(i_pcplot).Cmax    = sscanf(Fields{4}, '%f');
      Cdata.PcolorPlots(i_pcplot).Cticks  = sscanf(regexprep(Fields{5}, ':', ' '), '%f');
      Cdata.PcolorPlots(i_pcplot).Cscale  = Fields{6};
      Cdata.PcolorPlots(i_pcplot).Title   = ParseTitle(Fields{7});
      Cdata.PcolorPlots(i_pcplot).Xlabel  = regexprep(Fields{8}, '_', ' ');
      Cdata.PcolorPlots(i_pcplot).Ylabel  = regexprep(Fields{9}, '_', ' ');
      Cdata.PcolorPlots(i_pcplot).OutFile = Fields{10};
    case 'ProfTsPlot:'
      i_pts_plot = i_pts_plot + 1;
      Cdata.ProfTsPlots(i_pts_plot).PSname      = Fields{2};
      Cdata.ProfTsPlots(i_pts_plot).Fprefix     = Fields{3};
      Cdata.ProfTsPlots(i_pts_plot).Var         = Fields{4};
      Cdata.ProfTsPlots(i_pts_plot).Title       = ParseTitle(Fields{5});
      Cdata.ProfTsPlots(i_pts_plot).Tlabel      = regexprep(Fields{6}, '_', ' ');
      Cdata.ProfTsPlots(i_pts_plot).Zlabel      = regexprep(Fields{7}, '_', ' ');
      Cdata.ProfTsPlots(i_pts_plot).Cmin        = sscanf(Fields{8}, '%f');
      Cdata.ProfTsPlots(i_pts_plot).Cmax        = sscanf(Fields{9}, '%f');
      Cdata.ProfTsPlots(i_pts_plot).Zmin        = sscanf(Fields{10}, '%f');
      Cdata.ProfTsPlots(i_pts_plot).Zmax        = sscanf(Fields{11}, '%f');
      Cdata.ProfTsPlots(i_pts_plot).Pspec       = Fields{12};
      Cdata.ProfTsPlots(i_pts_plot).Flevel      = sscanf(Fields{13}, '%i');
      Cdata.ProfTsPlots(i_pts_plot).Ptype       = Fields{14};
      Cdata.ProfTsPlots(i_pts_plot).OutFileBase = Fields{15};

      Cdata.ProfTsPlots(i_pts_plot).PSnum       = -1;
    case 'ProfRsPlot:'
      i_prs_plot = i_prs_plot + 1;
      Cdata.ProfRsPlots(i_prs_plot).PSname      = Fields{2};
      Cdata.ProfRsPlots(i_prs_plot).Fprefix     = Fields{3};
      Cdata.ProfRsPlots(i_prs_plot).Var         = Fields{4};
      Cdata.ProfRsPlots(i_prs_plot).Title       = ParseTitle(Fields{5});
      Cdata.ProfRsPlots(i_prs_plot).Tlabel      = regexprep(Fields{6}, '_', ' ');
      Cdata.ProfRsPlots(i_prs_plot).Zlabel      = regexprep(Fields{7}, '_', ' ');
      Cdata.ProfRsPlots(i_prs_plot).Cmin        = sscanf(Fields{8}, '%f');
      Cdata.ProfRsPlots(i_prs_plot).Cmax        = sscanf(Fields{9}, '%f');
      Cdata.ProfRsPlots(i_prs_plot).Zmin        = sscanf(Fields{10}, '%f');
      Cdata.ProfRsPlots(i_prs_plot).Zmax        = sscanf(Fields{11}, '%f');
      Cdata.ProfRsPlots(i_prs_plot).Pspec       = Fields{12};
      Cdata.ProfRsPlots(i_prs_plot).Flevel      = sscanf(Fields{13}, '%i');
      Cdata.ProfRsPlots(i_prs_plot).Ptype       = Fields{14};
      Cdata.ProfRsPlots(i_prs_plot).OutFileBase = Fields{15};

      Cdata.ProfRsPlots(i_prs_plot).PSnum       = -1;
    case 'AzavgEofPlot:'
      i_aeofplot = i_aeofplot + 1;
      Cdata.AzavgEofPlots(i_aeofplot).Fprefix = Fields{2};
      Cdata.AzavgEofPlots(i_aeofplot).Var     = Fields{3};
      Cdata.AzavgEofPlots(i_aeofplot).Name    = regexprep(Fields{4}, '_', ' ');
      Cdata.AzavgEofPlots(i_aeofplot).Units   = regexprep(Fields{5}, '_', ' ');
      Cdata.AzavgEofPlots(i_aeofplot).Num     = sscanf(Fields{6},  '%d');
      Cdata.AzavgEofPlots(i_aeofplot).Clim    = sscanf(Fields{7},  '%f');
      Cdata.AzavgEofPlots(i_aeofplot).Cinc    = sscanf(Fields{8},  '%f');
      Cdata.AzavgEofPlots(i_aeofplot).Rmin    = sscanf(Fields{9},  '%f');
      Cdata.AzavgEofPlots(i_aeofplot).Rmax    = sscanf(Fields{10}, '%f');
      Cdata.AzavgEofPlots(i_aeofplot).Zmin    = sscanf(Fields{11},  '%f');
      Cdata.AzavgEofPlots(i_aeofplot).Zmax    = sscanf(Fields{12}, '%f');
  end
end

% Make the association between TsavgPlots and the PlotSets
if (isfield(Cdata, 'TsavgPlots'))
  Cdata.TsavgPlots = AssociateStructs( Cdata.TsavgPlots, Cdata.PlotSets, 'PS', 'PlotSet', 'TsavgPlot' ); 
end

% Make the association between DistPlots and the PlotSets
if (isfield(Cdata, 'DistPlots'))
  Cdata.DistPlots = AssociateStructs( Cdata.DistPlots, Cdata.PlotSets, 'PS', 'PlotSet', 'DistPlot' ); 
end
  
% Make the association between LinePlots and the PlotSets, PlotVars, PlotDselects
if (isfield(Cdata, 'LinePlots'))
  Cdata.LinePlots = AssociateStructs( Cdata.LinePlots, Cdata.PlotSets, 'PS', 'PlotSet', 'LinePlot' ); 
  Cdata.LinePlots = AssociateStructs( Cdata.LinePlots, Cdata.PlotVars, 'XV', 'PlotVar', 'LinePlot' ); 
  Cdata.LinePlots = AssociateStructs( Cdata.LinePlots, Cdata.PlotVars, 'YV', 'PlotVar', 'LinePlot' ); 
  Cdata.LinePlots = AssociateStructs( Cdata.LinePlots, Cdata.PlotDselects, 'DS', 'PlotDselect', 'LinePlot' ); 
end
  
% Make the association between ProfPlots and the PlotSets
if (isfield(Cdata, 'ProfPlots'))
  Cdata.ProfPlots = AssociateStructs( Cdata.ProfPlots, Cdata.PlotSets, 'PS', 'PlotSet', 'ProfPlot' ); 
end

% Make the association between ProfTsPlots and the PlotSets
if (isfield(Cdata, 'ProfTsPlots'))
  Cdata.ProfTsPlots = AssociateStructs( Cdata.ProfTsPlots, Cdata.PlotSets, 'PS', 'PlotSet', 'ProfTsPlot' ); 
end

% Make the association between ProfRsPlots and the PlotSets
if (isfield(Cdata, 'ProfRsPlots'))
  Cdata.ProfRsPlots = AssociateStructs( Cdata.ProfRsPlots, Cdata.PlotSets, 'PS', 'PlotSet', 'ProfRsPlot' ); 
end

% Make the association between HmeasSlicePlots and the PlotSets
if (isfield(Cdata, 'HmeasSlicePlots'))
  Cdata.HmeasSlicePlots = AssociateStructs( Cdata.HmeasSlicePlots, Cdata.PlotSets, 'PS', 'PlotSet', 'HmeasSlicePlot' ); 
end

% Make the association between Slope1dPlots and the PlotSets
if (isfield(Cdata, 'Slope1dPlots'))
  Cdata.Slope1dPlots = AssociateStructs( Cdata.Slope1dPlots, Cdata.PlotSets, 'PS', 'PlotSet', 'Slope1dPlot' ); 
end

% Make the association between Smeas and the PlotSets
if (isfield(Cdata, 'Smeas'))
  Cdata.Smeas = AssociateStructs( Cdata.Smeas, Cdata.PlotSets, 'PS', 'PlotSet', 'Smeas' ); 
end

% Make the association between Smeas1d and the PlotSets
if (isfield(Cdata, 'Smeas1d'))
  Cdata.Smeas1d = AssociateStructs( Cdata.Smeas1d, Cdata.PlotSets, 'PS', 'PlotSet', 'Smeas1d' ); 
end

% Make the association between CompCounts and the PlotSets
if (isfield(Cdata, 'CompCounts'))
  Cdata.CompCounts = AssociateStructs( Cdata.CompCounts, Cdata.PlotSets, 'PS', 'PlotSet', 'CompCounts' ); 
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AssociateStructs
%
% This function will take names specified in one struct array (Base) and find
% these names in another struct array (List), then record the matching indeces
% in the output structure OutStruct. Base is first copied to OutStruct so that
% the caller can replace the struct passed in by List. This gets around the fact
% that MATLAB uses call-by-value making it so that you cannot directly change
% an argument to this function.
%
% Type is used to enable Base to have multiple associations in it. The
% naming scheme for the structur elements needs to adhere to the following
% convection for this routine to work.
%
% 1) Base is an array of structures with an element named '<Type>name'
%     eg, if Type is 'PS', then Base need an element named 'PSname'
% 2) The value stored in Base(i).<Type>name will be searched for in List
%    where List is an array of structures with and element named 'Name'
%     eg, List(:).Name will be checked to see if it matches Base(i).PSname
% 3) The index in List where the match occurred will be entered into an element
%    in OutStruct (copy of Base) named '<Type>num'. A zero is assigned if
%    no match was made.
%     eg, if List(3).Name matched Base(2).PSname, then OutStruct(2).PSnum will
%     be set to 3.
%
% A copy of the input array with the new field (containing the index of the match)
% is returned to facilitate the replacement of Base in the calling routine.
% Since MATLAB uses pass-by-value it does not work to do the assignment
% here in this routine.
% 
function [ OutStruct ] = AssociateStructs(Base, List, Type, Ltype, Btype)

  % copy Base so the caller can replace Base
  OutStruct = Base;

  for i = 1:length(Base)
    BaseName = sprintf('%sname', Type);
    BaseNum  = sprintf('%snum', Type);

    TestName = Base(i).(BaseName);
    TestArray = { List(:).Name };
    
    Index = find(strcmp(TestName, TestArray));
    if (isempty(Index))
      fprintf('WARNING: Could not find a match for %s "%s" specified in %s number %d\n', Ltype, TestName, Btype, i);
      Index = 0;
    end

    OutStruct(i).(BaseNum) = Index;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ParseTitle
%
% This function will read the string given for a title and
% split it up into a main title and a set of panel markers
% (if given). The panel markers are denoted by:
%
%   PANEL:a:b:c:d
%
% where each string between ':'s is a successive label for a panel
% in the overall figure.
function [ Title ] = ParseTitle(Tstr)
  %
  % In general, underscores are replaced with spaces. Once
  % this is done, if the first token starts with "PANEL:", then
  % this token is further broken down into a set of one character
  % labels for each panel.

  Fields = regexp(Tstr, '_', 'split');
  if (regexp(Fields{1}, '^PANEL:'))
    % Split the first string into a list of panel markers
    Pmkrs = regexp(Fields{1}, ':', 'split');

    % At this point we want to toss the first element of Fields and Pmkrs. The
    % first element of Fields is 'PANEL:.*', and the first element of Pmkrs is 'PANEL'.
    if (length(Fields) >= 2)
      Title.Main = Fields{2};
      for i = 3:length(Fields)
        Title.Main = sprintf('%s %s', Title.Main, Fields{i});
      end
    else
      Title.Main = '';
    end
    for i = 2:length(Pmkrs)
      Title.Pmarkers{i-1} = Pmkrs{i};
    end
  else
    % no panel markers, use title as is
    Title.Main = Tstr;
    Title.Pmarkers = {}; % empty cell array
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Line2Fields
%
% This function will parse the string given in Line using the
% delimiter character given in Delim and separate Line into fields
% which get passed back to the caller in the array Fields.
%
function [ Fields ] = Line2Fields ( Line, Delim )

Remain = Line;

i = 1;
while (~isempty(Remain))
    [ Fields{i}, Remain ] = strtok(Remain, Delim);
    i = i + 1;
end

% Fields will end up with blank string in its last entry as
% a result of the above loop. If there is more than one element
% in Fields, then strip off this last entry since it is an artifact
if ((strcmp(Fields{end},'')) && (length(Fields) > 1))
  Fields = Fields(1:end-1);
end

end
