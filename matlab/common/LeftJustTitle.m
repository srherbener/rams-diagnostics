function [ ] = LeftJustTitle( Title )
%LeftJustTitle left justify title
%
%  Title is a handle to a title object
%

  % The title is in a box that adjusts to the amount of characters in
  % the title. Ie, it doesn't do any good to do Left/Center/Right alignment.
  % But, the entire box can be moved to the left side of the plot.

  set(Title, 'Units', 'Normalized');
  set(Title, 'HorizontalAlignment', 'Left');
  Tpos = get(Title, 'Position');
  Tpos(1) = 0; % line up with left edge of plot area
  set(Title, 'Position', Tpos);

end
