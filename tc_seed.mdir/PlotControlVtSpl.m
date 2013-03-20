function [ ] = PlotControlVtSpl(ConfigFile)
% PlotControlVtSpl function to plot max sfc Vt and SPL of control simulation
%

[ Config ] = ReadConfig(ConfigFile);

UndefVal = Config.UndefVal;
SpinUpCase = Config.SpinUpCase;
ControlCase = Config.ControlCase;
AzavgDir = Config.AzavgDir;
PlotDir = Config.PlotDir;
SstVal = Config.SstVal;

% Make sure PlotDir exists
if (exist(PlotDir, 'dir') ~= 7)
    mkdir(PlotDir);
end

% Need to stitch together the first 23 hours of the spin up run with
% the remainder (24 through 144 hours) of the control run in order to
% show the complete control run.
%
% VT will be (r,z,t), SPL will be (p,t)

% read in the Vt and time coordinate data
Hfile = sprintf('%s/speed_t_%s.h5', AzavgDir, SpinUpCase);
VT_SPINUP = squeeze(hdf5read(Hfile, '/speed_t'));
T = squeeze(hdf5read(Hfile, '/t_coords'));
Hfile = sprintf('%s/speed_t_%s.h5', AzavgDir, ControlCase);
VT_CNTL = squeeze(hdf5read(Hfile, '/speed_t'));
VT = cat(3, VT_SPINUP(:,:,1:24), VT_CNTL);

% read in the SPL
Hfile = sprintf('%s/sea_press_%s.h5', AzavgDir, SpinUpCase);
SPL_SPINUP = squeeze(hdf5read(Hfile, '/sea_press'));
Hfile = sprintf('%s/sea_press_%s.h5', AzavgDir, ControlCase);
SPL_CNTL = squeeze(hdf5read(Hfile, '/sea_press'));
SPL = cat(2, SPL_SPINUP(:,1:24), SPL_CNTL);

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
% Times = zeros(1,Ntsteps);
% it = 0;
% for i = 1:Ntsteps
%   % Times will be in hours
%   Times(i) = Tstart + ((i-1)*Tinc);
% 
%   if (mod(i-TsStart,TsPeriod) == 0)
%     it = it + 1;
%     % returns day and time strings
%     [ Ds, Ts ] = TimeToString(BT(1), BT(2), BT(3), BT(4), BT(5), BT(6), Times(i));
% 
%     Tticks(it) = Times(i);
%     Tlabels{it} = sprintf('%s/%s', Ds, Ts);
%   end
% end

% Use simulation time
Times = T / 3600;

% smooth out the time series
Flen = 5;
MAX_VT  = SmoothFillTseries(MAX_VT, Ntsteps, Flen);
MIN_SPL = SmoothFillTseries(MIN_SPL, Ntsteps, Flen);

% Plot
Lwidth = 2;
Fsize = 20;
Pfile = sprintf('%s/ControlVtSpl.fig', PlotDir);
Ptitle = sprintf('Control Storm Evolution');
%Xlabel = sprintf('Local Time, Starting at %s', StartTime);
Xlabel = sprintf('Simulation Time (hr)');
Y1label = sprintf('Maximum surface Vt (m/s)');
Y2label = sprintf('Minimum sea level pressure (mb)');

AxisPos = [ 0.11 0.15 0.75 0.75 ] ;
Xlims = [ 0 150 ];

% Set up plot controls according to SST
switch SstVal
  case 26
    Y1lims = [ 0 50 ];
    Y2lims = [ 970 1020 ];
    Tbefore = 24;
    Tduring = 54;
    Tafter = 84;
  case 27
    Y1lims = [ 0 60 ];
    Y2lims = [ 940 1020 ];
    Tbefore = 22;
    Tduring = 40;
    Tafter = 78;
  case 28
    Y1lims = [ 0 80 ];
    Y2lims = [ 900 1020 ];
    Tbefore = 24;
    Tduring = 42;
    Tafter = 60;
  otherwise
    % SST == 26 setup
    Y1lims = [ 0 50 ];
    Y2lims = [ 970 1020 ];
    Tbefore = 24;
    Tduring = 54;
    Tafter = 84;
end

% locations on Vt that correpsond to times relative to the
% rapid instensification phase
%   before - 24hrs
%   during - 54hrs
%   after  - 84hrs
%
% each index in Vt is one hour and the first index is zero so
% need to add one to the time to get the corresponding index.

Xbefore = [ Tbefore Tbefore ];
TxtBefore = sprintf('\\leftarrow %d hr', Tbefore);

Xduring = [ Tduring Tduring ];
TxtDuring = sprintf('\\leftarrow %d hr', Tduring);

Xafter = [ Tafter Tafter ];
TxtAfter = sprintf('\\leftarrow %d hr', Tafter);

Fig = figure;

% data
[ AX, H1, H2 ] = plotyy(Times, MAX_VT, Times, MIN_SPL);
set(H1, 'LineWidth', Lwidth, 'Color', 'k', 'LineStyle', '-');
set(H2, 'LineWidth', Lwidth, 'Color', 'k', 'LineStyle', '--');

% axes
axes(AX(1));
set(gca, 'FontSize', Fsize);
%set(gca, 'XTick', Tticks);
%set(gca, 'XTickLabel', Tlabels);
set(gca, 'Position', AxisPos);
set(gca, 'YColor', 'k');
xlim(Xlims);
ylabel(Y1label);
ylim(Y1lims);

axes(AX(2));
set(gca, 'FontSize', Fsize);
%set(gca, 'XTick', Tticks);
%set(gca, 'XTickLabel', Tlabels);
set(gca, 'Position', AxisPos);
set(gca, 'YColor', 'k');
xlim(Xlims);
ylabel(Y2label);
ylim(Y2lims);

%title(Ptitle);
xlabel(Xlabel);

% markers
Ytext = Y2lims(1) + ((Y2lims(2) - Y2lims(1))*0.95);
TfontSize = 14;
line(Xbefore,  Y2lims, 'LineStyle', ':', 'Color', 'k', 'LineWidth', Lwidth/2);
text(Xbefore(1), Ytext, TxtBefore, 'FontSize', TfontSize);
%line(Xduring, Y2lims, 'LineStyle', ':', 'Color', 'k', 'LineWidth', Lwidth/2);
%text(Xduring(1), Ytext, TxtDuring, 'FontSize', TfontSize);
%line(Xafter,  Y2lims, 'LineStyle', ':', 'Color', 'k', 'LineWidth', Lwidth/2);
%text(Xafter(1), Ytext, TxtAfter, 'FontSize', TfontSize);

saveas(Fig, Pfile);
close(Fig);
end
