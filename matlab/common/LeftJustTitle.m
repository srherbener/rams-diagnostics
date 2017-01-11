function [ ] = LeftJustTitle( Paxes, Title )
%LeftJustTitle left justify title
%
%  Paxes is a handle to the axes object of which contains the title
%  Title is a handle to the title object
%

  % The title is in a box that adjusts to the amount of characters in
  % the title. Ie, it doesn't do any good to do Left/Center/Right alignment.
  % But, the entire box can be moved to the left side of the plot.
  %
  % Using normalized units would seem to make sense since you would always
  % shift the title box to the x position of zero. This works in the matlab
  % figure but doesn't always get rendered correctly in a picture
  % file (jpg, png, etc.). It seems much more reliable to keep the units
  % of the title position as "data", and use the xlim property of the axes
  % to find the minimum (left side) x value.

  Xlim = get(Paxes, 'Xlim');

  set(Title, 'HorizontalAlignment', 'Left');
  Tpos = get(Title, 'Position');
  Tpos(1) = Xlim(1); % line up with left edge of plot area
  set(Title, 'Position', Tpos);

end
