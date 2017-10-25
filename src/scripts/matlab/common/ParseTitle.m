%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ParseTitle
%
% This function will read the string given for a title and
% split it up into a main title and a set of panel markers
% (if given). The panel markers are denoted by:
%
%   PANEL:a:b:c:d
%
% where each string between ':'s is a successive label for a panel
% in the overall figure.
function [ Title ] = ParseTitle(Tstr)
  %
  % In general, underscores are replaced with spaces. Once
  % this is done, if the first token starts with "PANEL:", then
  % this token is further broken down into a set of one character
  % labels for each panel.

  if (strcmp(Tstr, '@'))
      Tstr = '';
      Fields{1} = '';
  else
      Fields = regexp(Tstr, '@', 'split');
  end
  if (regexp(Fields{1}, '^PANEL:'))
    % Split the first string into a list of panel markers
    Pmkrs = regexp(Fields{1}, ':', 'split');

    % At this point we want to toss the first element of Fields and Pmkrs. The
    % first element of Fields is 'PANEL:.*', and the first element of Pmkrs is 'PANEL'.
    if (length(Fields) >= 2)
      Title.Main = Fields{2};
      for i = 3:length(Fields)
        Title.Main = sprintf('%s %s', Title.Main, Fields{i});
      end
    else
      Title.Main = '';
    end
    for i = 2:length(Pmkrs)
      Title.Pmarkers{i-1} = Pmkrs{i};
    end
  else
    % no panel markers, use title as is
    Title.Main = regexprep(Tstr, '@', ' ');
    Title.Pmarkers = {}; % empty cell array
  end
end
