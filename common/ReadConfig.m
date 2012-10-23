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
i_aeofplot = 0;
i_tsplot = 0;
i_dplot = 0;
i_2dplot = 0;
i_hmp3d = 0;
i_pplot = 0;
i_pset = 0;
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
    case 'EofDir:'
      Cdata.EofDir = Fields{2};
    case 'ControlCase:'
      Cdata.ControlCase = Fields{2};
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
      Cdata.AzavgEof(i_aeof).Name = Fields{2};
      Cdata.AzavgEof(i_aeof).Rvar = Fields{3};
      Cdata.AzavgEof(i_aeof).Rmin = sscanf(Fields{4}, '%f');
      Cdata.AzavgEof(i_aeof).Rmax = sscanf(Fields{5}, '%f');
      Cdata.AzavgEof(i_aeof).Zmin = sscanf(Fields{6}, '%f');
      Cdata.AzavgEof(i_aeof).Zmax = sscanf(Fields{7}, '%f');
      Cdata.AzavgEof(i_aeof).Tmin = sscanf(Fields{8}, '%f');
      Cdata.AzavgEof(i_aeof).Tmax = sscanf(Fields{9}, '%f');
    case 'Hmeas:'
      i_hmeas = i_hmeas + 1;
      Cdata.Hmeas(i_hmeas).Name    = Fields{2};
      Cdata.Hmeas(i_hmeas).InDir   = Fields{3};
      Cdata.Hmeas(i_hmeas).Fprefix = Fields{4};
      Cdata.Hmeas(i_hmeas).Rvar    = Fields{5};
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
        Cdata.PlotSets(i_pset).Cases(ips).Cname = Fields{j};
        Cdata.PlotSets(i_pset).Cases(ips).Legend = regexprep(Fields{j+1}, '_', ' ');
        j = j + 2;
      end
    case 'TsavgPlot:'
      i_tsplot = i_tsplot + 1;
      Cdata.TsavgPlots(i_tsplot).PSname  = Fields{2};
      Cdata.TsavgPlots(i_tsplot).Var     = Fields{3};
      Cdata.TsavgPlots(i_tsplot).Type    = Fields{4};
      Cdata.TsavgPlots(i_tsplot).Name    = regexprep(Fields{5}, '_', ' ');
      Cdata.TsavgPlots(i_tsplot).Units   = regexprep(Fields{6}, '_', ' ');
      Cdata.TsavgPlots(i_tsplot).Title   = regexprep(Fields{7}, '_', ' ');
      Cdata.TsavgPlots(i_tsplot).LegLoc  = Fields{8};
      Cdata.TsavgPlots(i_tsplot).Ymin    = sscanf(Fields{9},  '%f');
      Cdata.TsavgPlots(i_tsplot).Ymax    = sscanf(Fields{10},  '%f');
      Cdata.TsavgPlots(i_tsplot).OutFile = Fields{11};
    case 'DistPlot:'
      i_dplot = i_dplot + 1;
      Cdata.DistPlots(i_dplot).PSname  = Fields{2};
      Cdata.DistPlots(i_dplot).Var     = Fields{3};
      Cdata.DistPlots(i_dplot).Title   = regexprep(Fields{4}, '_', ' ');
      Cdata.DistPlots(i_dplot).Xlabel  = regexprep(Fields{5}, '_', ' ');
      Cdata.DistPlots(i_dplot).Ylabel  = regexprep(Fields{6}, '_', ' ');
      Cdata.DistPlots(i_dplot).LegLoc  = Fields{7};
      Cdata.DistPlots(i_dplot).Ymin    = sscanf(Fields{8},  '%f');
      Cdata.DistPlots(i_dplot).Ymax    = sscanf(Fields{9},  '%f');
      Cdata.DistPlots(i_dplot).OutFile = Fields{10};
    case 'TwoDimPlot:'
      i_2dplot = i_2dplot + 1;
      Cdata.TwoDimPlots(i_2dplot).PSname  = Fields{2};
      Cdata.TwoDimPlots(i_2dplot).Xvar    = Fields{3};
      Cdata.TwoDimPlots(i_2dplot).Yvar    = Fields{4};
      Cdata.TwoDimPlots(i_2dplot).Title   = regexprep(Fields{5}, '_', ' ');
      Cdata.TwoDimPlots(i_2dplot).Xlabel  = regexprep(Fields{6}, '_', ' ');
      Cdata.TwoDimPlots(i_2dplot).Ylabel  = regexprep(Fields{7}, '_', ' ');
      Cdata.TwoDimPlots(i_2dplot).LegLoc  = Fields{8};
      Cdata.TwoDimPlots(i_2dplot).OutFile = Fields{9};
    case 'HmeasPlot3d:'
      i_hmp3d = i_hmp3d + 1;
      Cdata.HmeasPlot3d(i_hmp3d).Name    = Fields{2};
      Cdata.HmeasPlot3d(i_hmp3d).Var     = Fields{3};
      Cdata.HmeasPlot3d(i_hmp3d).Fprefix = Fields{4};
      Cdata.HmeasPlot3d(i_hmp3d).Descrip = regexprep(Fields{5}, '_', ' ');
      Cdata.HmeasPlot3d(i_hmp3d).Units   = Fields{6};
      Cdata.HmeasPlot3d(i_hmp3d).Vtype   = Fields{7};
      Cdata.HmeasPlot3d(i_hmp3d).Isurf   = sscanf(Fields{8}, '%f');
      Cdata.HmeasPlot3d(i_hmp3d).Rmin    = sscanf(Fields{9}, '%f');
      Cdata.HmeasPlot3d(i_hmp3d).Rmax    = sscanf(Fields{10},'%f');
      Cdata.HmeasPlot3d(i_hmp3d).Zmin    = sscanf(Fields{11},'%f');
      Cdata.HmeasPlot3d(i_hmp3d).Zmax    = sscanf(Fields{12},'%f');
      Cdata.HmeasPlot3d(i_hmp3d).Tmin    = sscanf(Fields{13},'%f');
      Cdata.HmeasPlot3d(i_hmp3d).Tmax    = sscanf(Fields{14},'%f');
      Cdata.HmeasPlot3d(i_hmp3d).Cmin    = sscanf(Fields{15},'%f');
      Cdata.HmeasPlot3d(i_hmp3d).Cmax    = sscanf(Fields{16},'%f');
    case 'ProfPlot:'
      i_pplot = i_pplot + 1;
      Cdata.ProfPlots(i_pplot).PSname  = Fields{2};
      Cdata.ProfPlots(i_pplot).Var     = Fields{3};
      Cdata.ProfPlots(i_pplot).Type    = Fields{4};
      Cdata.ProfPlots(i_pplot).Title   = regexprep(Fields{5}, '_', ' ');
      Cdata.ProfPlots(i_pplot).Xlabel  = regexprep(Fields{6}, '_', ' ');
      Cdata.ProfPlots(i_pplot).Zlabel  = regexprep(Fields{7}, '_', ' ');
      Cdata.ProfPlots(i_pplot).LegLoc  = Fields{8};
      Cdata.ProfPlots(i_pplot).Xspec   = sscanf(Fields{9}, '%f:%f:%f', 3);
      Cdata.ProfPlots(i_pplot).Zmin    = sscanf(Fields{10}, '%f');
      Cdata.ProfPlots(i_pplot).Zmax    = sscanf(Fields{11}, '%f');
      Cdata.ProfPlots(i_pplot).OutFile = Fields{12};
    case 'AzavgEofPlot:'
      i_aeofplot = i_aeofplot + 1;
      Cdata.AzavgEofPlots(i_aeofplot).Var   = Fields{2};
      Cdata.AzavgEofPlots(i_aeofplot).Name  = regexprep(Fields{3}, '_', ' ');
      Cdata.AzavgEofPlots(i_aeofplot).Units = regexprep(Fields{4}, '_', ' ');
      Cdata.AzavgEofPlots(i_aeofplot).Num   = sscanf(Fields{5},  '%d');
      Cdata.AzavgEofPlots(i_aeofplot).Clim  = sscanf(Fields{6},  '%f');
      Cdata.AzavgEofPlots(i_aeofplot).Cinc  = sscanf(Fields{7},  '%f');
      Cdata.AzavgEofPlots(i_aeofplot).Rmin  = sscanf(Fields{8},  '%f');
      Cdata.AzavgEofPlots(i_aeofplot).Rmax  = sscanf(Fields{9},  '%f');
      Cdata.AzavgEofPlots(i_aeofplot).Zmin  = sscanf(Fields{10},  '%f');
      Cdata.AzavgEofPlots(i_aeofplot).Zmax  = sscanf(Fields{11}, '%f');
  end
end

% Make the association between TsavgPlots and the PlotSets
if (isfield(Cdata, 'TsavgPlots'))
  for itp = 1:length(Cdata.TsavgPlots)
    PlotSetName = Cdata.TsavgPlots(itp).PSname;
    Match = 0;
    for ips = 1:length(Cdata.PlotSets)
      if (strcmp(Cdata.PlotSets(ips).Name, PlotSetName))
        Match = ips;
      end
    end
    Cdata.TsavgPlots(itp).PSnum = Match;
  
    if (Match == 0)
      fprintf('WARNING: Could not find a match for PlotSet "%s" specified in TsavgPlot number %d\n', Cdata.TsavgPlots(itp).PSname, itp);
    end
  end
end

% Make the association between DistPlots and the PlotSets
if (isfield(Cdata, 'DistPlots'))
  for itp = 1:length(Cdata.DistPlots)
    PlotSetName = Cdata.DistPlots(itp).PSname;
    Match = 0;
    for ips = 1:length(Cdata.PlotSets)
      if (strcmp(Cdata.PlotSets(ips).Name, PlotSetName))
        Match = ips;
      end
    end
    Cdata.DistPlots(itp).PSnum = Match;
  
    if (Match == 0)
      fprintf('WARNING: Could not find a match for PlotSet "%s" specified in DistPlot number %d\n', Cdata.DistPlots(itp).PSname, itp);
    end
  end
end
  
% Make the association between TwoDimPlots and the PlotSets
if (isfield(Cdata, 'TwoDimPlots'))
  for itp = 1:length(Cdata.TwoDimPlots)
    PlotSetName = Cdata.TwoDimPlots(itp).PSname;
    Match = 0;
    for ips = 1:length(Cdata.PlotSets)
      if (strcmp(Cdata.PlotSets(ips).Name, PlotSetName))
        Match = ips;
      end
    end
    Cdata.TwoDimPlots(itp).PSnum = Match;
  
    if (Match == 0)
      fprintf('WARNING: Could not find a match for PlotSet "%s" specified in TwoDimPlot number %d\n', Cdata.TwoDimPlots(itp).PSname, itp);
    end
  end
end
  
% Make the association between ProfPlots and the PlotSets
if (isfield(Cdata, 'ProfPlots'))
  for itp = 1:length(Cdata.ProfPlots)
    PlotSetName = Cdata.ProfPlots(itp).PSname;
    Match = 0;
    for ips = 1:length(Cdata.PlotSets)
      if (strcmp(Cdata.PlotSets(ips).Name, PlotSetName))
        Match = ips;
      end
    end
    Cdata.ProfPlots(itp).PSnum = Match;
  
    if (Match == 0)
      fprintf('WARNING: Could not find a match for PlotSet "%s" specified in ProfPlot number %d\n', Cdata.ProfPlots(itp).PSname, itp);
    end
  end
end

end
