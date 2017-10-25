function [ ] = PlotInitialConditions(ConfigFile)
% PlotInitialConditions Plot the RAMS initialization profiles of Theta, Qt, U and V

Config = ReadConfig(ConfigFile);

Pdir = Config.PlotDir;
if (exist(Pdir, 'dir') ~= 7)
  mkdir(Pdir);
end

% From RAMS output
%--------REFERENCE STATE at I,J=(   0,   0)   SFC ELEV (M)=    0.0-------------
%
%    Z      U01D    V01D    DN01D    PI01D    PRESS    TH01D    THD      RT01D
%   (m)     (m/s)   (m/s)  (kg/m3)  (J/kgK)    (Pa)     (K)     (K)     (kg/kg)
REF_STATE = [
   -50.0  -10.85   -1.97   1.186   1009.93  102080.1  298.07  295.75   0.01288
    50.0  -10.85   -1.97   1.177   1006.64  100922.3  298.07  295.75   0.01288
   150.0  -10.55   -1.90   1.167   1003.35   99773.8  298.03  295.75   0.01263
   250.0  -10.25   -1.75   1.158   1000.06   98634.6  298.02  295.75   0.01261
   350.0   -9.95   -1.61   1.148    996.77   97504.6  298.02  295.75   0.01258
   450.0   -9.65   -1.46   1.139    993.49   96384.0  298.02  295.75   0.01256
   550.0   -9.35   -1.32   1.129    990.20   95272.5  298.01  295.75   0.01254
   650.0   -9.05   -1.17   1.120    986.91   94170.2  298.01  295.75   0.01251
   750.0   -8.75   -1.00   1.110    983.62   93077.5  298.20  296.12   0.01150
   850.0   -8.45   -0.87   1.100    980.34   91994.5  298.42  296.37   0.01131
   950.0   -8.15   -0.74   1.090    977.05   90921.4  298.64  296.62   0.01112
  1050.0   -7.85   -0.60   1.080    973.77   89857.9  298.85  296.87   0.01092
  1150.0   -7.55   -0.47   1.070    970.49   88804.2  299.07  297.12   0.01073
  1250.0   -7.25   -0.34   1.061    967.22   87760.0  299.29  297.37   0.01054
  1350.0   -6.95   -0.21   1.051    963.95   86725.5  299.50  297.62   0.01035
  1450.0   -6.59   -0.08   1.036    960.68   85702.9  301.21  299.55   0.00910
  1550.0   -6.17    0.05   1.017    957.45   84697.1  304.41  303.15   0.00680
  1650.0   -5.75    0.18   0.998    954.24   83710.2  307.59  306.75   0.00450
  1750.0   -5.46    0.29   0.988    951.06   82737.1  307.94  307.10   0.00450
  1850.0   -5.18    0.40   0.979    947.88   81773.2  308.29  307.45   0.00450
  1950.0   -4.89    0.51   0.970    944.70   80818.4  308.64  307.80   0.00450
  2050.0   -4.60    0.62   0.960    941.53   79872.6  309.00  308.15   0.00450
  2150.0   -4.31    0.73   0.951    938.36   78935.9  309.35  308.50   0.00450
  2250.0   -4.03    0.84   0.942    935.19   78008.1  309.70  308.85   0.00450
  2350.0   -3.74    0.95   0.933    932.03   77089.1  310.05  309.20   0.00450
  2450.0   -3.45    1.05   0.924    928.87   76179.0  310.40  309.55   0.00450
  2550.0   -3.16    1.16   0.915    925.72   75277.5  310.75  309.90   0.00450
  2650.0   -2.88    1.27   0.907    922.56   74384.7  311.10  310.25   0.00450
  2750.0   -2.59    1.38   0.898    919.42   73500.5  311.45  310.60   0.00450
  2850.0   -2.30    1.49   0.889    916.27   72624.8  311.80  310.95   0.00450
  2950.0   -2.02    1.60   0.881    913.13   71757.5  312.15  311.30   0.00450
  3050.0   -1.73    1.71   0.872    909.99   70898.6  312.51  311.65   0.00450
  3150.0   -1.44    1.82   0.864    906.86   70048.0  312.86  312.00   0.00450
  3250.0   -1.15    1.93   0.855    903.73   69205.7  313.21  312.35   0.00450
  3350.0   -0.87    2.04   0.847    900.60   68371.6  313.56  312.70   0.00450
  3450.0   -0.58    2.15   0.839    897.48   67545.6  313.91  313.05   0.00450
  3550.0   -0.29    2.26   0.831    894.36   66727.7  314.26  313.40   0.00450
  3650.0   -0.01    2.37   0.822    891.24   65917.8  314.61  313.75   0.00450
  3750.0    0.28    2.48   0.814    888.13   65115.8  314.96  314.10   0.00450
  3850.0    0.57    2.59   0.806    885.02   64321.6  315.31  314.45   0.00450
  3950.0    0.86    2.70   0.798    881.91   63535.3  315.66  314.80   0.00450
  ];

% Grab the profiles from the RAMS report
% Discard the first row from the table above since it is a RAMS level below the surface
Z   = squeeze(REF_STATE(2:end,1)) ./ 1000; % km
U   = squeeze(REF_STATE(2:end,2));         % m/s
V   = squeeze(REF_STATE(2:end,3));         % m/s
DN  = squeeze(REF_STATE(2:end,4));         % kg/m3
PI  = squeeze(REF_STATE(2:end,5));         % J/(kg K)
P   = squeeze(REF_STATE(2:end,6)) ./ 100;  % mb
TH  = squeeze(REF_STATE(2:end,7));         % K
THD = squeeze(REF_STATE(2:end,8));         % K
QT  = squeeze(REF_STATE(2:end,9)) .* 1000; % g/kg

% Create three panels
%   Theta
%   Qt
%   U and V

Fsize = 20;

OutFile = sprintf('%s/InitialConditions.jpg', Pdir);

Fig = figure;

% Put Theta on the left, label the y-axis with height
subplot(1,3,1);
plot(TH,Z,'LineWidth', 2)
xlim([ 295 320 ]);
set(gca, 'FontSize', Fsize);

LeftJustifyTitle(title('a)'));
xlabel('\theta (K)');
ylabel('Height (km)');

% Put Qt in the middle, shut off y-axis labeling
subplot(1,3,2);
plot(QT,Z,'LineWidth', 2)
xlim([ 0 15 ]);
set(gca, 'FontSize', Fsize);
set(gca, 'YtickLabel', {});

LeftJustifyTitle(title('b)'));
xlabel('q_t (g kg^-^1)');

% Put U and V on the right, shut off y-axis labeling
subplot(1,3,3);
plot([U V],Z,'LineWidth', 2)
xlim([ -15 5 ]);
set(gca, 'FontSize', Fsize);
set(gca, 'YtickLabel', {});

LeftJustifyTitle(title('c)'));
xlabel('U, V (m s^-^1)');

legend({ 'U' 'V' }, 'FontSize', 14, 'Location', 'NorthWest');

saveas(Fig, OutFile);
close(Fig);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to line up the title with the left side of the plot area
function [] = LeftJustifyTitle(T)
  set(T, 'Units', 'Normalized');
  set(T, 'HorizontalAlignment', 'Left');
  Tpos = get(T, 'Position');
  Tpos(1) = 0;  % line up with left edge of plot area
  set(T,'Position', Tpos);
end
