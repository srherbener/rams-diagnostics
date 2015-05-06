function [ Cdata ] = ReadFigureConfig ( Cfile )
% ReadFigureConfig read configuration file for constructing figures

  % Parse the input config file into statements (lists of tokens)
  Cstmts = ReadPrepConfigFile(Cfile);

  % Read through the statements and build the output structure containing
  % the description of figures to create.
  i_pset = 0;
  i_pdat = 0;
  i_paxes = 0;
  i_fpan = 0;

  % Provide nan as a defalut for the undefined value
  Cdata.UndefVal = nan;
  
  for is = 1:length(Cstmts)
    Fields = Cstmts{is};
    
    switch(Fields{1})

      case 'UndefVal:'
        Cdata.UndefVal = sscanf(Fields{2}, '%f');
        
      case 'PlotSet:'
        i_pset = i_pset + 1;
        Cdata.PlotSets(i_pset).Name   = Fields{2};
        Cdata.PlotSets(i_pset).Nlines = sscanf(Fields{3}, '%d');
        j = 4; % next field
        for ipl = 1:Cdata.PlotSets(i_pset).Nlines
          Cdata.PlotSets(i_pset).Lines(ipl).PDname  = Fields{j};
          Cdata.PlotSets(i_pset).Lines(ipl).PDnum   = -1;
          Cdata.PlotSets(i_pset).Lines(ipl).Legend  = regexprep(Fields{j+1}, '@', ' ');
          Cdata.PlotSets(i_pset).Lines(ipl).Lwidth  = sscanf(Fields{j+2}, '%f');
          Cdata.PlotSets(i_pset).Lines(ipl).Lcolor  = Fields{j+3};
          Cdata.PlotSets(i_pset).Lines(ipl).Lstyle  = Fields{j+4};
          Cdata.PlotSets(i_pset).Lines(ipl).Lgscale = sscanf(Fields{j+5}, '%f');
          Cdata.PlotSets(i_pset).Lines(ipl).Xzoom   = sscanf(Fields{j+6}, '%f');
          Cdata.PlotSets(i_pset).Lines(ipl).Yzoom   = sscanf(Fields{j+7}, '%f');
          j = j + 8;
        end
        
      case 'PlotData:'
        i_pdat = i_pdat + 1;
        Cdata.PlotData(i_pdat).Name   = Fields{2};
        % X data
        Cdata.PlotData(i_pdat).Xvar    = Fields{3};
        Cdata.PlotData(i_pdat).Xfile   = Fields{4};
        Cdata.PlotData(i_pdat).XSelect = regexprep(Fields{5}, '@', ' ');;
        Cdata.PlotData(i_pdat).Xscale  = sscanf(Fields{6}, '%f');
        Cdata.PlotData(i_pdat).Xoffset = sscanf(Fields{7}, '%f');
        % Y data
        Cdata.PlotData(i_pdat).Yvar    = Fields{8};
        Cdata.PlotData(i_pdat).Yfile   = Fields{9};
        Cdata.PlotData(i_pdat).YSelect = regexprep(Fields{10}, '@', ' ');;
        Cdata.PlotData(i_pdat).Yscale  = sscanf(Fields{11},  '%f');
        Cdata.PlotData(i_pdat).Yoffset = sscanf(Fields{12}, '%f');
        
      case 'PlotAxes:'
        i_paxes = i_paxes + 1;
        Cdata.PlotAxes(i_paxes).Name    = Fields{2};
        Cdata.PlotAxes(i_paxes).Fsize   = sscanf(Fields{3}, '%f');
        Cdata.PlotAxes(i_paxes).Lwidth  = sscanf(Fields{4}, '%f');
        Cdata.PlotAxes(i_paxes).Tlength = eval(regexprep(Fields{5}, '_', ' '));
        % X axis
        Cdata.PlotAxes(i_paxes).Xlabel   = regexprep(Fields{6}, '@', ' ');
        Cdata.PlotAxes(i_paxes).Xunits   = regexprep(Fields{7}, '@', ' ');
        Cdata.PlotAxes(i_paxes).Xmin     = sscanf(Fields{8}, '%f');
        Cdata.PlotAxes(i_paxes).Xmax     = sscanf(Fields{9}, '%f');
        Cdata.PlotAxes(i_paxes).Xscale   = Fields{10};
        Cdata.PlotAxes(i_paxes).Xticks   = eval(regexprep(Fields{11}, '_', ' '));
        % Y axis
        Cdata.PlotAxes(i_paxes).Ylabel   = regexprep(Fields{12}, '@', ' ');
        Cdata.PlotAxes(i_paxes).Yunits   = regexprep(Fields{13}, '@', ' ');
        Cdata.PlotAxes(i_paxes).Ymin     = sscanf(Fields{14}, '%f');
        Cdata.PlotAxes(i_paxes).Ymax     = sscanf(Fields{15}, '%f');
        Cdata.PlotAxes(i_paxes).Yscale   = Fields{16};
        Cdata.PlotAxes(i_paxes).Yticks   = eval(regexprep(Fields{17}, '_', ' '));
        
      case 'FigPanel:'
        i_fpan = i_fpan + 1;
        Cdata.FigPanels(i_fpan).Name     = Fields{2};
        Cdata.FigPanels(i_fpan).PSname   = Fields{3};
        Cdata.FigPanels(i_fpan).PSnum    = -1;
        Cdata.FigPanels(i_fpan).PAname   = Fields{4};
        Cdata.FigPanels(i_fpan).PAnum    = -1;
        Cdata.FigPanels(i_fpan).XAshow   = sscanf(Fields{5}, '%d');
        Cdata.FigPanels(i_fpan).YAshow   = sscanf(Fields{6}, '%d');
        Cdata.FigPanels(i_fpan).Smooth   = Fields{7};
        Cdata.FigPanels(i_fpan).Title    = ParseTitle(Fields{8});
        Cdata.FigPanels(i_fpan).LegLoc   = Fields{9};
        Cdata.FigPanels(i_fpan).LegFsize = sscanf(Fields{10}, '%d');
    end
  end

  % Make the association between PlotSets and PlotData
  if (isfield(Cdata, 'PlotSets'))
    for i_pset = 1:length(Cdata.PlotSets)
      Cdata.PlotSets(i_pset).Lines = AssociateStructs( Cdata.PlotSets(i_pset).Lines, Cdata.PlotData, 'PD', 'PlotData', 'PlotSets' ); 
    end
  end

  % Make the association between FigPanels and the PlotSets, PlotData andPlotAxes
  if (isfield(Cdata, 'FigPanels'))
    Cdata.FigPanels = AssociateStructs( Cdata.FigPanels, Cdata.PlotSets,     'PS', 'PlotSet',     'LinePlot' ); 
    Cdata.FigPanels = AssociateStructs( Cdata.FigPanels, Cdata.PlotAxes,     'PA', 'PlotAxes',    'LinePlot' ); 
  end
end
