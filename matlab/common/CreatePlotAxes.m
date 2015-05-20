function [ Axes ] = CreatePlotAxes( Prows, Pcols, Ploc, AxesSpecs )
%CreatePlotAxes create axes for a plot
%   Axes - handle to a figure that was opened by the caller
%
%   AxesProps - structure containing a list of axis property names and
%               associated values that are desired to be set
%

  % select the passed in figure
  Axes = subplot(Prows, Pcols, Ploc);

  AxesProps = AxesSpecs.Props;
  Nprops = length(AxesProps);
  for i = 1:Nprops
    set(Axes, AxesProps(i).Name, AxesProps(i).Val);
  end

  % Title and axes labels
  T = title(AxesSpecs.Title);
  xlabel(AxesSpecs.Xlabel);
  ylabel(AxesSpecs.Ylabel);

  % If have a panel marker, left justify the title
  if (AxesSpecs.Panel)
    LeftJustTitle(T);
  end

end
