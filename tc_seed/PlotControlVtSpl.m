function [ ] = PlotControlVtSpl(ConfigFile)
% PlotControlVtSpl function to plot max sfc Vt and SPL of control simulation
%

[ Config ] = ReadConfig(ConfigFile);

UndefVal = Config.UndefVal;
ControlCase = Config.ControlCase;
AzavgDir = Config.AzavgDir;
PlotDir = Config.PlotDir;

% read in the Vt and time coordinate data
Hfile = sprintf('%s/speed_t_%s.h5', AzavgDir, ControlCase);
VT = squeeze(hdf5read(Hfile, '/speed_t'));
T = squeeze(hdf5read(Hfile, '/t_coords'));

% read in the SPL
Hfile = sprintf('%s/sea_press_%s.h5', AzavgDir, ControlCase);
SPL = squeeze(hdf5read(Hfile, '/sea_press'));

% VT will be (r,z,t), SPL will be (p,t)
% Create a time series of the max sfc Vt and one of the min SPL
% z == 1 is below sfc, so use z == 2 for sfc Vt
VT = squeeze(VT(:,2,:)); % take the z == 2 level

MIN_SPL = min(SPL,[],1);
MAX_VT  = max(VT,[],1);

% Figure out the times and time labels
BT(1) = 1998;
BT(2) = 8;
BT(3) = 21;
BT(4) = 21;
BT(5) = 0;
BT(6) = 0;
BaseTime = datenum(BT(1), BT(2), BT(3), BT(4), BT(5), BT(6));
StartTime = datestr(BaseTime, 'mm/dd/yyyy HH:MM');

Ntsteps  = length(T);
Tstart   = 0;
Tinc     = 1;
TsStart  = 10;
TsPeriod = 30;

% Generate the time values
Times = zeros(1,Ntsteps);
it = 0;
for i = 1:Ntsteps
  % Times will be in hours
  Times(i) = Tstart + ((i-1)*Tinc);

  if (mod(i-TsStart,TsPeriod) == 0)
    it = it + 1;
    % returns day and time strings
    [ Ds, Ts ] = TimeToString(BT(1), BT(2), BT(3), BT(4), BT(5), BT(6), Times(i));

    Tticks(it) = Times(i);
    Tlabels{it} = sprintf('%s/%s', Ds, Ts);
  end
end

% smooth out the time series
Flen = 5;
MAX_VT  = SmoothFillTseries(MAX_VT, Ntsteps, Flen);
MIN_SPL = SmoothFillTseries(MIN_SPL, Ntsteps, Flen);

% Plot
Lwidth = 2;
Fsize = 20;
Pfile = sprintf('%s/CntlVtSpl.jpg', PlotDir);
Ptitle = sprintf('Control Run Storm Development');
Xlabel = sprintf('Local Time, Starting at %s', StartTime);
Y1label = sprintf('Maximum surface Vt (m/s)');
Y2label = sprintf('Minimum sea level pressure (mb)');
AxisPos = [ 0.11 0.15 0.75 0.75 ] ;

Fig = figure;

% data
[ AX, H1, H2 ] = plotyy(Times, MAX_VT, Times, MIN_SPL);
set(H1, 'LineWidth', Lwidth);
set(H2, 'LineWidth', Lwidth);

% axes
axes(AX(1));
set(gca, 'FontSize', Fsize);
set(gca, 'XTick', Tticks);
set(gca, 'XTickLabel', Tlabels);
set(gca, 'Position', AxisPos);
ylabel(Y1label);

axes(AX(2));
set(gca, 'FontSize', Fsize);
set(gca, 'XTick', Tticks);
set(gca, 'XTickLabel', Tlabels);
set(gca, 'Position', AxisPos);
ylabel(Y2label);

title(Ptitle);
xlabel(Xlabel);

saveas(Fig, Pfile);
close(Fig);
end
