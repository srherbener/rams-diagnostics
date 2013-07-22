function screen2jpeg(fig, filename)
%SCREEN2JPEG Generate a JPEG file of the given figure with
% dimensions consistent with the figure's screen dimensions.
%
% SCREEN2JPEG(fig, 'filename') saves the given figure to the
% JPEG file "filename".
%
% Sean P. McCarthy
% Copyright (c) 1984-98 by MathWorks, Inc. All Rights Reserved
figure(fig);

if nargin < 1
error('Not enough input arguments!')
end

oldscreenunits = get(gcf,'Units');
oldpaperunits = get(gcf,'PaperUnits');
oldpaperpos = get(gcf,'PaperPosition');
set(gcf,'Units','pixels');
scrpos = get(gcf,'Position');
newpos = scrpos/100;
set(gcf,'PaperUnits','inches',...
'PaperPosition',newpos)
print('-djpeg', filename, '-r100');
drawnow
set(gcf,'Units',oldscreenunits,...
'PaperUnits',oldpaperunits,...
'PaperPosition',oldpaperpos);

end