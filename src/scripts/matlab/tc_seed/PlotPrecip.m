function [ ] = PlotPrecip()
% PlotPrecip function to plot precip rate
%

  Pdir = 'Plots';
  Ddir = 'DIAGS';
  
  CaseList = {
    { 'TCS_SD_C0100' 'C0100' 'blue'   2  '-' }
    { 'TCS_SD_C0500' 'C0500' 'green'  2  '-' }
    { 'TCS_SD_C1000' 'C1000' 'red'    2  '-' }
    { 'TCS_SD_C2000' 'C2000' 'black'  2  '-' }
    };
  
  Nc = length(CaseList);
  
  Hvar = '/pcprate';
  Rvar = '/radius';

  % Set up axes for the plot
  AxesSpecs.Title = '';
  AxesSpecs.Xlabel = 'Radius (km)';
  AxesSpecs.Ylabel = 'Precip Rate (mm h^-^1)';

  i_ap = 1;

  % set up axes properties
  AxesSpecs.Props(i_ap).Name = 'FontSize';
  AxesSpecs.Props(i_ap).Val  = 25;
  i_ap = i_ap + 1;

  % Read in line data for the plot
  for icase = 1:Nc
    Case           = CaseList{icase}{1};
    LegText{icase} = CaseList{icase}{2};

    PlotData.Lcolor{icase}  = CaseList{icase}{3};
    PlotData.Lwidth{icase}  = CaseList{icase}{4};
    PlotData.Lstyle{icase}  = CaseList{icase}{5};
  
    % read in the Vt and coordinate data
    % TE data is (x,y,z,t) where x is radius
    Hfile = sprintf('%s/pcprate_%s_SS.h5', Ddir, Case);
    fprintf('Reading: %s (%s, %s)\n', Hfile, Hvar, Rvar);
    PR = squeeze(h5read(Hfile, Hvar)); % mm/h
    R  = squeeze(h5read(Hfile, Rvar)); % km
  
    % Build data for the plot
    PlotData.X{icase} = R;
    PlotData.Y{icase} = PR;
  end
  fprintf('\n');

  % create plot and save in jpeg file
  OutFile = sprintf('%s/TimeAvgPrecipRate_SD.jpg', Pdir);
  fprintf('Writing: %s\n', OutFile);

  Fig = figure;
  MakeLinePlot(Fig, AxesSpecs, PlotData);
  legend(LegText, 'Location', 'NorthEast', 'FontSize', 15);

  saveas(Fig, OutFile);
  close(Fig);

  fprintf('\n');

end
