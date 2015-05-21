function [ Axes ] = SetPlotAxes( Axes, AxesSpecs )
%SetPlotAxes set axes properties
%   Axes - handle to a figure that was opened by the caller
%
%   AxesProps - structure containing a list of axis property names and
%               associated values that are desired to be set
%

  % select the passed in axes
  axes(Axes);

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

  % If have color axis limits, set them
  if (~isempty(AxesSpecs.Clims))
    caxis(AxesSpecs.Clims);
  end

end
