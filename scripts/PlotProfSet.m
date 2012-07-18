function [ ] = PlotProfSet( Xvals, Zvals, Profs, Xlabel, Zlabel, Ptitle, Lcolors, LegText, LegLoc, OutFile )
%PlotProfSet Plot a set of vertical profiles on the same panel
%   This function will take data contained in Profs and plot them on a
%   single panel line plot. Each row of Profs is a separate profile.
%
%   OutFile is the path to the file that contains the image of the plot.

Fig = figure;

% Need to establish an axis style before calling "hold on" (for the
% subsequent plots). Plot the first profile outside the loop, issue a "hold
% on" and then plot the remainder profiles.

% Want to rotate plot 90 deg counter clockwise which turns out to be +90
% deg around the Z axis. After rotation, Xvals will apply to the x-axis and
% Zvals will apply to the y-axis.

[ Nprofs, Npts ] = size(Profs);
Lwidth = 2;

P = plot(Zvals,Profs(1,:),'LineWidth',Lwidth,'Color',char(Lcolors(1)));
rotate(P, [ 0 0 1 ], 90, [ 0 0 0 ]);
xlim([ min(Xvals) max(Xvals) ]);
ylim([ min(Zvals) max(Zvals) ]);
set (gca, 'FontSize', 20);
hold on;

for i = 2:Nprofs
    P = plot(Zvals,Profs(i,:),'LineWidth',Lwidth,'Color',char(Lcolors(i)));
    rotate(P, [ 0 0 1 ], 90, [ 0 0 0 ]);
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

