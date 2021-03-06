function [ ] = PlotProfSet( Xvals, Zvals, Profs, Xlabel, Zlabel, Ptitle, Pmarkers, Lcolors, Lstyles, Lgscales, LegText, LegLoc, Ptype, OutFile )
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
  Lwidth = 4;
  Fsize = 35;
  LegFsize = 20;
end

PanelTitle = ~isempty(Pmarkers);
if (PanelTitle)
  Ptitle = sprintf('(%s) %s', Pmarkers{1}, Ptitle);
else
  Lwidth = 2;
  Fsize = 30;
  LegFsize = 20;
end

if (strcmp(Lcolors{1}, 'k'))
  Lcolor = [ 1 1 1 ] * Lgscales(1);
else
  Lcolor = str2rgb(Lcolors{1});
end
plot(Profs(1,:),Zvals,'Color', Lcolor, 'LineStyle', Lstyles{1},'LineWidth',Lwidth);
xlim([ min(Xvals) max(Xvals) ]);
ylim([ min(Zvals) max(Zvals) ]);
set (gca, 'FontSize', Fsize);
set(gca, 'LineWidth', 2);
set(gca, 'TickLength', [ 0.025 0.025 ]);
hold on;

for i = 2:Nprofs
     if (strcmp(Lcolors{i}, 'k'))
       Lcolor = [ 1 1 1 ] * Lgscales(i);
     else
       Lcolor = str2rgb(Lcolors{i});
     end
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

if (PanelTitle)
  % Fix up the positioning
  Ppos = get(gca, 'Position'); % position of plot area
  Ppos(1) = Ppos(1) * 1.00;
  Ppos(2) = Ppos(2) * 0.95;
  Ppos(3) = Ppos(3) * 0.90;
  Ppos(4) = Ppos(4) * 0.90;
  set(gca, 'Position', Ppos);
end

saveas(Fig, OutFile);

hold off;
close(Fig);

end

