function [ ] = PlotProfSet( Xvals, Zvals, Profs, Xlabel, Zlabel, Ptitle, Lstyles, Lgscales, LegText, LegLoc, Ptype, OutFile )
%PlotProfSet Plot a set of vertical profiles on the same panel
%   This function will take data contained in Profs and plot them on a
%   single panel line plot. Each row of Profs is a separate profile.
%
%   OutFile is the path to the file that contains the image of the plot.

Fig = figure;

% Need to establish an axis style before calling "hold on" (for the
% subsequent plots). Plot the first profile outside the loop, issue a "hold
% on" and then plot the remainder profiles.

[ Nprofs, Npts ] = size(Profs);
if (strcmp(Ptype, 'test'))
  Lwidth = 1;
  Fsize = 16;
  LegFsize = 12;
else
  Lwidth = 3;
  Fsize = 45;
  LegFsize = 25;
end

Yticks = [ 2 6 10 14 ];
Yticklabels = { '2' '6' '10' '14' };

PanelTitle = false;
if (regexp(Ptitle, '^PANEL:'))
    Ptitle = regexprep(Ptitle, '^PANEL:', '');
    if (regexp(Ptitle, ' '))
        Ptitle = regexprep(Ptitle, ' ', ') ');
    else
        Ptitle = sprintf('%s)', Ptitle);
    end
    PanelTitle = true;
end

Lcolor = [ 1 1 1 ] * Lgscales(1);
plot(Profs(1,:),Zvals,'Color', Lcolor, 'LineStyle', Lstyles{1},'LineWidth',Lwidth);
xlim([ min(Xvals) max(Xvals) ]);
ylim([ min(Zvals) max(Zvals) ]);
set (gca, 'FontSize', Fsize);
set(gca, 'YTick', Yticks);
set(gca, 'YTickLabel', Yticklabels);
hold on;

for i = 2:Nprofs
     Lcolor = [ 1 1 1 ] * Lgscales(i);
     plot(Profs(i,:),Zvals,'Color', Lcolor, 'LineStyle', Lstyles{i},'LineWidth',Lwidth);
end

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
else
    title(Ptitle);
end
xlabel(Xlabel);
ylabel(Zlabel);
legend(LegText, 'Location', LegLoc, 'FontSize', LegFsize);
legend boxoff;

% Fix up the positioning
Ppos = get(gca, 'Position'); % position of plot area
Ppos(1) = Ppos(1) * 1.00;
Ppos(2) = Ppos(2) * 1.00;
Ppos(3) = Ppos(3) * 0.85;
Ppos(4) = Ppos(4) * 0.85;
set(gca, 'Position', Ppos);

saveas(Fig, OutFile);

hold off;
close(Fig);

end

