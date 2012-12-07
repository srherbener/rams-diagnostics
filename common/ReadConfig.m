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
i_aeofplot = 0;
i_tsplot = 0;
i_dplot = 0;
i_lplot = 0;
i_hmp3d = 0;
i_hmslice = 0;
i_pplot = 0;
i_pts_plot = 0;
i_pset = 0;
i_pvar = 0;
i_sspec = 0;
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
        Cdata.PlotSets(i_pset).Cases(ips).Lspec = Fields{j+2};
        j = j + 3;
      end
    case 'PlotVar:'
      i_pvar = i_pvar + 1;
      Cdata.PlotVars(i_pvar).Name    = Fields{2};
      Cdata.PlotVars(i_pvar).Var     = Fields{3};
      Cdata.PlotVars(i_pvar).Fprefix = Fields{4};
      Cdata.PlotVars(i_pvar).Label   = regexprep(Fields{5}, '_', ' ');;
      Cdata.PlotVars(i_pvar).Units   = regexprep(Fields{6}, '_', ' ');;
      Cdata.PlotVars(i_pvar).Min     = sscanf(Fields{7}, '%f');
      Cdata.PlotVars(i_pvar).Max     = sscanf(Fields{8}, '%f');
      Cdata.PlotVars(i_pvar).Scale   = sscanf(Fields{9}, '%f');
    case 'Sspec:'
      i_sspec = i_sspec + 1;
      Cdata.Sspecs(i_sspec).Name = Fields{2};
      Cdata.Sspecs(i_sspec).Xmin = sscanf(Fields{3}, '%f');
      Cdata.Sspecs(i_sspec).Xmax = sscanf(Fields{4}, '%f');
      Cdata.Sspecs(i_sspec).Ymin = sscanf(Fields{5}, '%f');
      Cdata.Sspecs(i_sspec).Ymax = sscanf(Fields{6}, '%f');
      Cdata.Sspecs(i_sspec).Zmin = sscanf(Fields{7}, '%f');
      Cdata.Sspecs(i_sspec).Zmax = sscanf(Fields{8}, '%f');
      Cdata.Sspecs(i_sspec).Tmin = sscanf(Fields{9}, '%f');
      Cdata.Sspecs(i_sspec).Tmax = sscanf(Fields{10}, '%f');
    case 'TsavgPlot:'
      i_tsplot = i_tsplot + 1;
      Cdata.TsavgPlots(i_tsplot).PSname  = Fields{2};
      Cdata.TsavgPlots(i_tsplot).Fprefix = Fields{3};
      Cdata.TsavgPlots(i_tsplot).Var     = Fields{4};
      Cdata.TsavgPlots(i_tsplot).Ptype   = Fields{5};
      Cdata.TsavgPlots(i_tsplot).Ttype   = Fields{6};
      Cdata.TsavgPlots(i_tsplot).Name    = regexprep(Fields{7}, '_', ' ');
      Cdata.TsavgPlots(i_tsplot).Units   = regexprep(Fields{8}, '_', ' ');
      Cdata.TsavgPlots(i_tsplot).Title   = regexprep(Fields{9}, '_', ' ');
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
      Cdata.DistPlots(i_dplot).Title   = regexprep(Fields{4}, '_', ' ');
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
      Cdata.LinePlots(i_lplot).XSname  = Fields{3};
      Cdata.LinePlots(i_lplot).YSname  = Fields{4};
      Cdata.LinePlots(i_lplot).SSname  = Fields{5};
      Cdata.LinePlots(i_lplot).Smooth  = Fields{6};
      Cdata.LinePlots(i_lplot).Title   = regexprep(Fields{7}, '_', ' ');
      Cdata.LinePlots(i_lplot).LegLoc  = Fields{8};
      Cdata.LinePlots(i_lplot).OutFile = Fields{9};

      Cdata.LinePlots(i_lplot).PSnum   = -1;
      Cdata.LinePlots(i_lplot).XSnum   = -1;
      Cdata.LinePlots(i_lplot).YSnum   = -1;
      Cdata.LinePlots(i_lplot).SSnum   = -1;
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
      Cdata.ProfPlots(i_pplot).Title   = regexprep(Fields{6}, '_', ' ');
      Cdata.ProfPlots(i_pplot).Xlabel  = regexprep(Fields{7}, '_', ' ');
      Cdata.ProfPlots(i_pplot).Zlabel  = regexprep(Fields{8}, '_', ' ');
      Cdata.ProfPlots(i_pplot).LegLoc  = Fields{9};
      Cdata.ProfPlots(i_pplot).Xspec   = sscanf(Fields{10}, '%f:%f:%f', 3);
      Cdata.ProfPlots(i_pplot).Zmin    = sscanf(Fields{11}, '%f');
      Cdata.ProfPlots(i_pplot).Zmax    = sscanf(Fields{12}, '%f');
      Cdata.ProfPlots(i_pplot).OutFile = Fields{13};

      Cdata.ProfPlots(i_pplot).PSnum   = -1;
    case 'ProfTsPlot:'
      i_pts_plot = i_pts_plot + 1;
      Cdata.ProfTsPlots(i_pts_plot).PSname      = Fields{2};
      Cdata.ProfTsPlots(i_pts_plot).Fprefix     = Fields{3};
      Cdata.ProfTsPlots(i_pts_plot).Var         = Fields{4};
      Cdata.ProfTsPlots(i_pts_plot).Title       = regexprep(Fields{5}, '_', ' ');
      Cdata.ProfTsPlots(i_pts_plot).Tlabel      = regexprep(Fields{6}, '_', ' ');
      Cdata.ProfTsPlots(i_pts_plot).Zlabel      = regexprep(Fields{7}, '_', ' ');
      Cdata.ProfTsPlots(i_pts_plot).Cmin        = sscanf(Fields{8}, '%f');
      Cdata.ProfTsPlots(i_pts_plot).Cmax        = sscanf(Fields{9}, '%f');
      Cdata.ProfTsPlots(i_pts_plot).Zmin        = sscanf(Fields{10}, '%f');
      Cdata.ProfTsPlots(i_pts_plot).Zmax        = sscanf(Fields{11}, '%f');
      Cdata.ProfTsPlots(i_pts_plot).Pspec       = Fields{12};
      Cdata.ProfTsPlots(i_pts_plot).OutFileBase = Fields{13};

      Cdata.ProfTsPlots(i_pts_plot).PSnum       = -1;
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
  
% Make the association between LinePlots and the PlotSets, PlotVars, Sspecs
if (isfield(Cdata, 'LinePlots'))
  Cdata.LinePlots = AssociateStructs( Cdata.LinePlots, Cdata.PlotSets, 'PS', 'PlotSet', 'LinePlot' ); 
  Cdata.LinePlots = AssociateStructs( Cdata.LinePlots, Cdata.PlotVars, 'XS', 'PlotVar', 'LinePlot' ); 
  Cdata.LinePlots = AssociateStructs( Cdata.LinePlots, Cdata.PlotVars, 'YS', 'PlotVar', 'LinePlot' ); 
  Cdata.LinePlots = AssociateStructs( Cdata.LinePlots, Cdata.Sspecs, 'SS', 'Sspec', 'LinePlot' ); 
end
  
% Make the association between ProfPlots and the PlotSets
if (isfield(Cdata, 'ProfPlots'))
  Cdata.ProfPlots = AssociateStructs( Cdata.ProfPlots, Cdata.PlotSets, 'PS', 'PlotSet', 'ProfPlot' ); 
end

% Make the association between ProfTsPlots and the PlotSets
if (isfield(Cdata, 'ProfTsPlots'))
  Cdata.ProfTsPlots = AssociateStructs( Cdata.ProfTsPlots, Cdata.PlotSets, 'PS', 'PlotSet', 'ProfTsPlot' ); 
end

% Make the association between HmeasSlicePlots and the PlotSets
if (isfield(Cdata, 'HmeasSlicePlots'))
  Cdata.HmeasSlicePlots = AssociateStructs( Cdata.HmeasSlicePlots, Cdata.PlotSets, 'PS', 'PlotSet', 'HmeasSlicePlot' ); 
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
    if (length(Index) == 0)
      fprintf('WARNING: Could not find a match for %s "%s" specified in %s number %d\n', Ltype, TestName, Btype, i);
      Index = 0;
    end

    OutStruct(i).(BaseNum) = Index;
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

end
