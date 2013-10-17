function [ ] = PlotControlVtSpl(ConfigFile, InType)
% PlotControlVtSpl function to plot max sfc Vt and SPL of control simulation
%

UseVt = strcmp(InType, 'vt');

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
if (UseVt)
  Hfile = sprintf('%s/speed_t_%s.h5', AzavgDir, SpinUpCase);
  VT_SPINUP = squeeze(hdf5read(Hfile, '/speed_t'));
else
  Hfile = sprintf('%s/speed10m_%s.h5', AzavgDir, SpinUpCase);
  VT_SPINUP = squeeze(hdf5read(Hfile, '/speed10m'));
end
T = squeeze(hdf5read(Hfile, '/t_coords'));
if (UseVt)
  Hfile = sprintf('%s/speed_t_%s.h5', AzavgDir, ControlCase);
  VT_CNTL = squeeze(hdf5read(Hfile, '/speed_t'));
else
  Hfile = sprintf('%s/speed10m_%s.h5', AzavgDir, ControlCase);
  VT_CNTL = squeeze(hdf5read(Hfile, '/speed10m'));
end

if (UseVt)
  VT = cat(3, VT_SPINUP(:,:,1:24), VT_CNTL);
else
  VT = cat(2, VT_SPINUP(:,1:24), VT_CNTL);
end

% read in the SPL
Hfile = sprintf('%s/sea_press_%s.h5', AzavgDir, SpinUpCase);
SPL_SPINUP = squeeze(hdf5read(Hfile, '/sea_press'));
Hfile = sprintf('%s/sea_press_%s.h5', AzavgDir, ControlCase);
SPL_CNTL = squeeze(hdf5read(Hfile, '/sea_press'));
SPL = cat(2, SPL_SPINUP(:,1:24), SPL_CNTL);

% Create a time series of the max sfc Vt and one of the min SPL
% z == 1 is below sfc, so use z == 2 for sfc Vt
if (UseVt)
  VT = squeeze(VT(:,2,:)); % take the z == 2 level
else
  VT = squeeze(VT(:,:));
end

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
Lwidth = 3;
Fsize = 35;
if (UseVt)
  Pfile = sprintf('%s/ControlVtSpl.jpg', PlotDir);
else
  Pfile = sprintf('%s/ControlVtSpl_s10.jpg', PlotDir);
end
Xlabel = sprintf('Simulation Time (hr)');
if (UseVt)
  Y1label = sprintf('Max Vt (m/s)');
else
  Y1label = sprintf('Max V_1_0_m (m/s)');
end
Y2label = sprintf('Min SLP (mb)');

AxisPos = [ 0.11 0.15 0.75 0.75 ] ;
Xlims = [ 0 150 ];

% Set up plot controls according to SST
switch SstVal
  case 26
    if (UseVt)
      Y1lims = [ 0 50 ];
      Y2lims = [ 970 1020 ];
    else
      Y1lims = [ 0 60 ];
      Y2lims = [ 960 1020 ];
    end
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
set(gca, 'LineWidth', 2);
set(gca, 'TickLength', [ 0.025 0.025 ]);
set(gca, 'Position', AxisPos);
set(gca, 'YColor', 'k');
set(gca, 'XTick', [ 50 100 ]);
if (UseVt)
  set(gca, 'YTick', [ 0 20 40 ]);
else
  set(gca, 'YTick', [ 0 20 40 60 ]);
end
xlim(Xlims);
ylabel(Y1label);
ylim(Y1lims);

axes(AX(2));
set(gca, 'FontSize', Fsize);
set(gca, 'LineWidth', 2);
set(gca, 'TickLength', [ 0.025 0.025 ]);
set(gca, 'Position', AxisPos);
set(gca, 'YColor', 'k');
set(gca, 'XTick', [ 50 100 ]);
if (UseVt)
  set(gca, 'YTick', [ 970 990 1010 ]);
else
  set(gca, 'YTick', [ 960 980 1000 1020 ]);
end
xlim(Xlims);
ylabel(Y2label);
ylim(Y2lims);

% The title is in a box that adjusts to the amount of characters in
% the title. Ie, it doesn't do any good to do Left/Center/Right
% alignment. But, the entire box can be moved to the left side of the
% plot.
T = title('(d)');
set(T, 'Units', 'Normalized');
set(T, 'HorizontalAlignment', 'Left');
Tpos = get(T, 'Position');
Tpos(1) = 0; % line up with left edge of plot area
set(T, 'Position', Tpos);
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

% fix position
Ppos = get(gca, 'Position');
Ppos(1) = Ppos(1) * 1.65;
Ppos(2) = Ppos(2) * 1.65;
Ppos(3) = Ppos(3) * 0.70;
Ppos(4) = Ppos(4) * 0.80;
set(gca, 'Position', Ppos);

saveas(Fig, Pfile);
close(Fig);
end
