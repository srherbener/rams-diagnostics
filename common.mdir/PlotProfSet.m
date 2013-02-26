function [ ] = PlotProfSet( Xvals, Zvals, Profs, Xlabel, Zlabel, Ptitle, Lspecs, LegText, LegLoc, OutFile )
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
Lwidth = 2;

plot(Profs(1,:),Zvals,Lspecs{1},'LineWidth',Lwidth);
xlim([ min(Xvals) max(Xvals) ]);
ylim([ min(Zvals) max(Zvals) ]);
set (gca, 'FontSize', 20);
hold on;

for i = 2:Nprofs
     plot(Profs(i,:),Zvals,Lspecs{i},'LineWidth',Lwidth);
end

title(Ptitle);
xlabel(Xlabel);
ylabel(Zlabel);
legend(LegText, 'Location', LegLoc);
legend boxoff;

saveas(Fig, OutFile);

hold off;
close(Fig);

end

