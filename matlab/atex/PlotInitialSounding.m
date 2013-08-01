function [ ] = PlotInitialSounding(ConfigFile)
% PlotInitialSounding Plot the ATEX sounding used in RAMS initialization

Config = ReadConfig(ConfigFile);

Pdir = Config.PlotDir;
if (exist(Pdir, 'dir') ~= 7)
  mkdir(Pdir);
end

% From RAMSIN
%     P           Z         T      Theta      U        V       Rt         RH
%     (Pa)       (m)       (K)      (K)      (m/s)    (m/s)   (kg/kg)     (%)
SOUNDING = [
    101500.0       0.0   297.01   295.75   -11.00    -2.00   0.01300     69.6
     99773.8     150.0   295.56   295.75   -10.55    -1.90   0.01263     72.8
     93622.6     700.0   290.23   295.75    -8.90    -1.10   0.01250     94.6
     93077.3     750.0   290.11   296.13    -8.75    -1.00   0.01150     87.2
     86211.6    1400.0   285.39   297.75    -6.80    -0.14   0.01025     98.1
     83709.0    1650.0   291.55   306.75    -5.75     0.18   0.00450     27.9
     63145.1    4000.0   276.19   314.98     1.00     2.75   0.00450     59.4
     ];

P     = squeeze(SOUNDING(:,1)) ./ 100;   % mb
Z     = squeeze(SOUNDING(:,2)) ./ 1000;  % km
T     = squeeze(SOUNDING(:,3)) - 273.15; % deg C
THETA = squeeze(SOUNDING(:,4));          % K
U     = squeeze(SOUNDING(:,5));          % m/s
V     = squeeze(SOUNDING(:,6));          % m/s
RT    = squeeze(SOUNDING(:,7));          % kg/kg
RH_PC = squeeze(SOUNDING(:,8));          % percent

RH = RH_PC ./ 100; % fraction
TD = T - ((100 - (RH_PC)) ./ 5);
WSPEED = sqrt(U.^2 + V.^2);
WDIR = atan2(V,U) .* (180/3.14159);

% WDIR is -180 to +180 using the trig convention of location on the axes
% (0 degrees is toward +x and increasing values go counter clockwise (CCW))
% A compass has zero toward +y (north) and increasing values go clockwise (CW)
WDIR(WDIR < 0) = WDIR(WDIR < 0) + 360;         % 0 to 360, increasing values are going CCW
WDIR = WDIR - 90;                              % -90 to 270, rotate so 0 points north
WDIR = 360 - WDIR;                             % 90 to 450, make increasing values go CW
WDIR(WDIR >= 360) = WDIR(WDIR >= 360) - 360;   % 0 to 360

% Use the metpack code
% zoom in on plot since the simulation top is 4000 m AGL
Fsize = 30;
OutFile = sprintf('%s/InitialSounding.jpg', Pdir);

Fig = figure;
set(gca, 'FontSize', Fsize);
skew_sounding(Fig, P, T, TD, RH, WSPEED, WDIR);

xlim([ -10 30 ]);
ylim([ 500 1050 ]);

saveas(Fig, OutFile);
close(Fig);

end
