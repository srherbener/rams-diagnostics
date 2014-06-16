function [ ] = GenFigures(PlotDir, FigureDir)
% GenFigures generate figures for the presentation and paper

  % each entry is one figure:
  %  {
  %  nrows
  %  ncols
  %    {
  %    PanelList
  %    }
  %  OutFile
  %  }
  %
  % PanelList is:
  %  { Panel 'Figure File' ShowX ShowY PanelMarker }
  %  ...
  %
  %  Panel is the selection in the subplot
  %  If figure file is 'legend', then creating a legend for the other panels
  %  ShowX, ShowY == 1 -> show the x axis and x label, otherwise hide these
  %
  FigSets = {
    {
    3
    3
      {
        { 1 'InvHeightStall_CO_S293_G10M5_ALL.fig' 0 1 'a' }
        { 2 'InvHeightStall_CO_S298_G10M5_ALL.fig' 0 0 'b' }
        { 3 'InvHeightStall_CO_S303_G10M5_ALL.fig' 0 0 'c' }
        { 4 'CloudFrac_CO_S293_G10M5_ALL.fig'      1 1 'd' }
        { 5 'CloudFrac_CO_S298_G10M5_ALL.fig'      1 0 'e' }
        { 6 'CloudFrac_CO_S303_G10M5_ALL.fig'      1 0 'f' }
      }
    'FigureTsBlStats_ALL.jpg'
    }

    };
  Nset = length(FigSets);

  fprintf('Generating figures:\n');

  for iset = 1:Nset
    Nrows     = FigSets{iset}{1};
    Ncols     = FigSets{iset}{2};
    PanelList = FigSets{iset}{3};
    OutFile   = sprintf('%s/%s', FigureDir, FigSets{iset}{4});

    % read in each plot file and place in panels of a subplot
    OutFig = figure('Visible', 'off');
    Nplots = length(PanelList);
    for iplot = 1:Nplots
      PanelNum = PanelList{iplot}{1};
      Pfile   = sprintf('%s/%s', PlotDir, PanelList{iplot}{2});
      ShowX   = PanelList{iplot}{3};
      ShowY   = PanelList{iplot}{4};
      Pmarker = PanelList{iplot}{5};

      fprintf('  Reading: %s\n', Pfile); 

      % input plot becomes a panel in the output subplot
      InFig  = openfig(Pfile, 'new', 'invisible');
      InAxes = gca;

      figure(OutFig);
      set(OutFig, 'Visible', 'off');
      OutAxes = subplot(Nrows, Ncols, PanelNum);

      CopyPlot(InAxes, OutAxes, ShowX, ShowY, Pmarker); 

      close(InFig);
    end
    fprintf('\n');


    fprintf('  Writing: %s\n', OutFile);
    saveas(OutFig, OutFile);
    fprintf('\n');

    close(OutFig);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ ] = CopyPlot(InAxes, OutAxes, ShowX, ShowY, Pmarker) 
% CopyPlot copy the graphics objects in InAxes to OutAxes

  % Copy the font size
  set(OutAxes, 'FontSize', get(InAxes, 'FontSize'));

  % Copy the lines
  Olist = findobj(InAxes, 'Type', 'line');
  copyobj(Olist, OutAxes);

  % Set the xlim and ylim
  set(OutAxes, 'Xlim', get(InAxes, 'Xlim'));
  set(OutAxes, 'Ylim', get(InAxes, 'Ylim'));

  % Set the visibility of the axes
  if (ShowX == 1)
    xlabel(get(get(InAxes, 'Xlabel'), 'String'));
  else
    set(OutAxes, 'Xtick', []);
  end
  if (ShowY == 1)
    ylabel(get(get(InAxes, 'Ylabel'), 'String'));
  else
    set(OutAxes, 'Ytick', []);
  end

  % Place the title
  PanelTitle = ~isempty(Pmarker);

  Ptitle = get(get(InAxes, 'Title'), 'String');
  if (PanelTitle)
    Ptitle = sprintf('(%s) %s', Pmarker, Ptitle);
  end

  if (~isempty(Ptitle))
    if (PanelTitle)
      % The title is in a box that adjusts to the amount of characters in
      % the title. Ie, it doesn't do any good to do Left/Center/Right
      % alignment. But, the entire box can be moved to the left side of the
      % plot.
      T = title(Ptitle);
      set(T, 'Units', 'Normalized');
      set(T, 'HorizontalAlignment', 'Left');
      Tpos = get(T, 'Position');
      Tpos(1) = 0; % line up with left edge of plot area
      set(T, 'Position', Tpos);
fprintf('Adjusting title\n');
    else
      title(Ptitle);
fprintf('Not adjusting title\n');
    end
  end

end
