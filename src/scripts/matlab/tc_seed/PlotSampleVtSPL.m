function [ ] = PlotSampleVtSPL(ConfigFile)
% PlotSampleVtSPL function to plot max sfc Vt and SPL of sample simulation
%

[ Config ] = ReadConfig(ConfigFile);

AzavgDir = Config.AzavgDir;
PlotDir = Config.PlotDir;

%DoingSD = false; % sea salt enabled
DoingSD = true; % sea salt disabled

if (DoingSD)
  Case1 = 'TCS_SD_C0100';
  Case2 = 'TCS_SD_C2000';
else
  Case1 = 'TCS_GN_C0100';
  Case2 = 'TCS_GN_C2000';
end

% for selection
Tmin = 24; 
Tmax = 90;

% for plot time axis
Xmin = 0;
Xmax = 100;

% Make sure PlotDir exists
if (exist(PlotDir, 'dir') ~= 7)
    mkdir(PlotDir);
end

% read in the Vt and time coordinate data
% VT is (r,z,t), use z = 2 which is first model level above ground
Hfile = sprintf('%s/speed_t_%s.h5', AzavgDir, Case1);
T = squeeze(hdf5read(Hfile, '/t_coords')) / 3600; % hr

T1 = find(T >= Tmin, 1, 'first');
T2 = find(T <= Tmax, 1, 'last');
Times = T(T1:T2);

VT = squeeze(hdf5read(Hfile, '/speed_t'));
VT = squeeze(VT(:,2,T1:T2)); % take the z = 2 level
MAX_VT1  = max(VT,[],1);

Hfile = sprintf('%s/speed_t_%s.h5', AzavgDir, Case2);
VT = squeeze(hdf5read(Hfile, '/speed_t'));
VT = squeeze(VT(:,2,T1:T2)); % take the z = 2 level
MAX_VT2  = max(VT,[],1);

% read in the SPL
% SPL is (p,t)
Hfile = sprintf('%s/sea_press_%s.h5', AzavgDir, Case1);
SPL = squeeze(hdf5read(Hfile, '/sea_press'));
SPL = squeeze(SPL(:,T1:T2));
MIN_SPL1 = min(SPL,[],1);

Hfile = sprintf('%s/sea_press_%s.h5', AzavgDir, Case2);
SPL = squeeze(hdf5read(Hfile, '/sea_press'));
SPL = squeeze(SPL(:,T1:T2));
MIN_SPL2 = min(SPL,[],1);


% smooth out the time series
Flen = 5;
Nt = length(Times);
MAX_VT1  = SmoothFillTseries(MAX_VT1, Nt, Flen);
MAX_VT2  = SmoothFillTseries(MAX_VT2, Nt, Flen);
MIN_SPL1 = SmoothFillTseries(MIN_SPL1, Nt, Flen);
MIN_SPL2 = SmoothFillTseries(MIN_SPL2, Nt, Flen);

% Plot
Lwidth = 2;
Fsize = 20;
if (DoingSD)
  Pfile = sprintf('%s/TCWS0513_SampleVtSpl.jpg', PlotDir);
  Ptitle = sprintf('Max Vt, Min SLP');
  Xlabel = sprintf('Simulation Time (hr)');
  Y1label = sprintf('Wind Speed (m/s)');
  Y2label = sprintf('SLP (mb)');
else
  Pfile = sprintf('%s/ams_SampleVtSpl.jpg', PlotDir);
  Ptitle = sprintf('Maximum Vt, Minimum Sea-level Pressure');
  Xlabel = sprintf('Simulation Time (hr)');
  Y1label = sprintf('Maximum surface Vt (m/s)');
  Y2label = sprintf('Minimum sea level pressure (mb)');
end

AxisPos = [ 0.11 0.15 0.75 0.75 ] ;
Xlims = [ Xmin Xmax ];
Y1lims = [ 0 60 ];
Y2lims = [ 960 1020 ];


% PLOT
Fig = figure;

% data
[ AX1, H1, H2 ] = plotyy(Times, MAX_VT1, Times, MIN_SPL1);
set(H1, 'LineWidth', Lwidth, 'Color', 'b', 'LineStyle', '-');
set(H2, 'LineWidth', Lwidth, 'Color', 'b', 'LineStyle', '--');

axes(AX1(1));
line(Times, MAX_VT2, 'Color', 'r', 'LineWidth', Lwidth, 'LineStyle', '-');
set(gca, 'FontSize', Fsize);
set(gca, 'Position', AxisPos);
set(gca, 'YColor', 'k');
xlim(Xlims);
ylabel(Y1label);
ylim(Y1lims);

axes(AX1(2));
line(Times, MIN_SPL2, 'Color', 'r', 'LineWidth', Lwidth, 'LineStyle', '--');
set(gca, 'FontSize', Fsize);
set(gca, 'Position', AxisPos);
set(gca, 'YColor', 'k');
xlim(Xlims);
ylabel(Y2label);
ylim(Y2lims);

title(Ptitle);
xlabel(Xlabel);

% In the SPL axes
text(Xmin+10, 967, 'CLEAN', 'FontSize', Fsize, 'Color', 'b');
text(Xmin+10, 963, 'POLLUTED', 'FontSize', Fsize, 'Color', 'r');

saveas(Fig, Pfile);
close(Fig);
end
