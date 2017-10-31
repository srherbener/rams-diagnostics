function [ ] = PlotPressWind(ConfigFile)


[ Config ] = ReadConfig(ConfigFile);

Pdir = Config.PlotDir;
if (exist(Pdir, 'dir') ~= 7)
    mkdir(Pdir);
end

% for the TS Debby simulations:
Times = [ 1:13 ];

NhcSLP = [ 1007 1007 1007 1007 1005 1003 1002 1001 1001 1000  999 1000 1000 ];
NhcWind  = [ 35   35   35   35   35   40   45   50   50   50   50   50   50   ]; 

SimSLP = [ 1006 1004 1004 1001 1002 1003 1002 1001 1000 1002 1001 1002 1001 ];
SimWind  = [ 36.2 41.4 47.4 44.7 55.3 55.7 54.1 53.2 48.3 44.1 44.7 45.4 45.2 ];

% plot
OutFile = sprintf('%s/TsDebbyPressWind.jpg', Pdir);

FigPressWind = figure;
subplot(2,1,1);
set(gca, 'FontSize', 18);
NhcSLP = plot(Times, NhcSLP, 'linewi', 1, 'color', 'b');
SimSLP = line(Times, SimSLP, 'linewi', 1, 'color', 'r');
title('TS Debby Minumum SLP');
%xlabel('Time');
set(gca,'xtick', [ 2, 6, 10 ]);
set(gca,'xticklabel', { 'Aug, 22, 00Z', 'Aug, 23, 00Z', 'Aug, 24, 00Z' });
ylabel('Pressure (mb)');
ylim([ 998 1008 ]);
legend([ NhcSLP SimSLP ], 'NHC Best Track', 'Simulated Track', 'Location', 'NorthEast');

subplot(2,1,2);
set(gca, 'FontSize', 18);
NhcWind = plot(Times, NhcWind, 'linewi', 1, 'color', 'b');
SimWind = line(Times, SimWind, 'linewi', 1, 'color', 'r');
title('TS Debby Maximum Surface Wind Speed');
xlabel('Time');
set(gca,'xtick', [ 2, 6, 10 ]);
set(gca,'xticklabel', { 'Aug, 22, 00Z', 'Aug, 23, 00Z', 'Aug, 24, 00Z' });
ylabel('Wind Speed (mph)');
ylim([ 30 60 ]);
legend([ NhcWind SimWind ], 'NHC Best Track', 'Simulated Track', 'Location', 'SouthEast');

fprintf('Writing: %s\n', OutFile);
saveas(FigPressWind, OutFile);
close(FigPressWind);
