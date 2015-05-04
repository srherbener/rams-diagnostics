function [ Cdata ] = ReadFigureConfig ( Cfile )
% ReadFigureConfig read configuration file for constructing figures

  % Parse the input config file into statements (lists of tokens)
  Cstmts = ReadPrepConfigFile(Cfile);

  % Read through the statements and build the output structure containing
  % the description of figures to create.
  i_pset = 0;
  i_pdat = 0;
  i_paxes = 0;
  i_pds = 0;
  i_lplot = 0;

  for is = 1:length(Cstmts)
    Fields = Cstmts{is};
    
    switch(Fields{1})
      case 'AzavgDir:'
        Cdata.AzavgDir = Fields{2};
        
      case 'TsavgDir:'
        Cdata.TsavgDir = Fields{2};
        
      case 'DiagDir:'
        Cdata.DiagDir = Fields{2};
        
      case 'PlotDir:'
        Cdata.PlotDir = Fields{2};
        
      case 'UndefVal:'
        Cdata.UndefVal = sscanf(Fields{2}, '%f');
        
      case 'PlotSet:'
        i_pset = i_pset + 1;
        Cdata.PlotSets(i_pset).Name   = Fields{2};
        Cdata.PlotSets(i_pset).Ncases = sscanf(Fields{3}, '%d');
        j = 4; % next field
        for ips = 1:Cdata.PlotSets(i_pset).Ncases
          Cdata.PlotSets(i_pset).Cases(ips).Cname   = Fields{j};
          Cdata.PlotSets(i_pset).Cases(ips).Legend  = regexprep(Fields{j+1}, '@', ' ');
          Cdata.PlotSets(i_pset).Cases(ips).Lcolor  = Fields{j+2};
          Cdata.PlotSets(i_pset).Cases(ips).Lstyle  = Fields{j+3};
          Cdata.PlotSets(i_pset).Cases(ips).Lgscale = sscanf(Fields{j+4}, '%f');
          Cdata.PlotSets(i_pset).Cases(ips).Xzoom   = sscanf(Fields{j+5}, '%f');
          Cdata.PlotSets(i_pset).Cases(ips).Yzoom   = sscanf(Fields{j+6}, '%f');
          j = j + 7;
        end
        
      case 'PlotData:'
        i_pdat = i_pdat + 1;
        Cdata.PlotData(i_pdat).Name   = Fields{2};
        % X data
        Cdata.PlotData(i_pdat).Xvar    = Fields{3};
        Cdata.PlotData(i_pdat).Xfile   = Fields{4};
        Cdata.PlotData(i_pdat).Xlabel  = regexprep(Fields{5}, '@', ' ');
        Cdata.PlotData(i_pdat).Xunits  = regexprep(Fields{6}, '@', ' ');
        Cdata.PlotData(i_pdat).Xscale  = sscanf(Fields{7}, '%f');
        Cdata.PlotData(i_pdat).Xoffset = sscanf(Fields{8}, '%f');
        % Y data
        Cdata.PlotData(i_pdat).Yvar    = Fields{9};
        Cdata.PlotData(i_pdat).Yfile   = Fields{10};
        Cdata.PlotData(i_pdat).Ylabel  = regexprep(Fields{11}, '@', ' ');
        Cdata.PlotData(i_pdat).Yunits  = regexprep(Fields{12}, '@', ' ');
        Cdata.PlotData(i_pdat).Yscale  = sscanf(Fields{13}, '%f');
        Cdata.PlotData(i_pdat).Yoffset = sscanf(Fields{14}, '%f');
        
      case 'PlotAxes:'
        i_paxes = i_paxes + 1;
        Cdata.PlotAxes(i_paxes).Name    = Fields{2};
        % X axis
        Cdata.PlotAxes(i_paxes).Xmin     = sscanf(Fields{3}, '%f');
        Cdata.PlotAxes(i_paxes).Xmax     = sscanf(Fields{4}, '%f');
        Cdata.PlotAxes(i_paxes).Xscale   = Fields{5};
        Cdata.PlotAxes(i_paxes).Xticks   = eval(regexprep(Fields{6}, '_', ' '));
        % Y axis
        Cdata.PlotAxes(i_paxes).Ymin     = sscanf(Fields{7}, '%f');
        Cdata.PlotAxes(i_paxes).Ymax     = sscanf(Fields{8}, '%f');
        Cdata.PlotAxes(i_paxes).Yscale   = Fields{9};
        Cdata.PlotAxes(i_paxes).Yticks   = eval(regexprep(Fields{10}, '_', ' '));
        
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
        
      case 'LinePlot:'
        i_lplot = i_lplot + 1;
        Cdata.LinePlots(i_lplot).PSname  = Fields{2};
        Cdata.LinePlots(i_lplot).PDname  = Fields{3};
        Cdata.LinePlots(i_lplot).PAname  = Fields{4};
        Cdata.LinePlots(i_lplot).XAshow  = sscanf(Fields{5}, '%d');
        Cdata.LinePlots(i_lplot).YAshow  = sscanf(Fields{6}, '%d');
        Cdata.LinePlots(i_lplot).DSname  = Fields{7};
        Cdata.LinePlots(i_lplot).Smooth  = Fields{8};
        Cdata.LinePlots(i_lplot).Title   = ParseTitle(Fields{9});
        Cdata.LinePlots(i_lplot).LegLoc  = Fields{10};
        Cdata.LinePlots(i_lplot).OutFile = Fields{11};

        Cdata.LinePlots(i_lplot).PSnum   = -1;
        Cdata.LinePlots(i_lplot).PDnum   = -1;
        Cdata.LinePlots(i_lplot).PAnum   = -1;
        Cdata.LinePlots(i_lplot).DSnum   = -1;
    end
  end
  
  % Make the association between LinePlots and the PlotSets, PlotData,
  % PlotAxes and PlotDselects
  if (isfield(Cdata, 'LinePlots'))
    Cdata.LinePlots = AssociateStructs( Cdata.LinePlots, Cdata.PlotSets,     'PS', 'PlotSet',     'LinePlot' ); 
    Cdata.LinePlots = AssociateStructs( Cdata.LinePlots, Cdata.PlotData,     'PD', 'PlotData',    'LinePlot' ); 
    Cdata.LinePlots = AssociateStructs( Cdata.LinePlots, Cdata.PlotAxes,     'PA', 'PlotAxes',    'LinePlot' ); 
    Cdata.LinePlots = AssociateStructs( Cdata.LinePlots, Cdata.PlotDselects, 'DS', 'PlotDselect', 'LinePlot' ); 
  end
end
