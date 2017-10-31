function [ ] = Sounding(Dummy)

% From RAMS log
%
%------------------------------SOUNDING INPUT----------------------------------
%
%       PS         HS       TS      THDS      US       VS       RTS     REL HUM
%      (Pa)       (m)      (K)      (K)      (m/s)    (m/s)   (kg/kg)     (%)
%

SOUNDING = [
    100000.0       0.0   300.00   300.00     0.00     0.00   0.01300     57.1
     97206.8     250.0   297.92   300.34     0.42     0.00   0.01300     63.0
     94473.4     500.0   295.96   300.81     0.83     0.00   0.01300     69.0
     91799.8     750.0   294.06   301.34     1.22     0.00   0.01300     75.5
     89184.7    1000.0   292.21   301.93     1.61     0.00   0.01211     76.9
     86626.8    1250.0   290.38   302.54     1.97     0.00   0.01095     75.8
     84125.2    1500.0   288.58   303.20     2.31     0.00   0.00989     74.8
     81679.4    1750.0   286.80   303.88     2.63     0.00   0.00893     73.7
     79289.0    2000.0   285.03   304.58     2.91     0.00   0.00806     72.5
     76953.1    2250.0   283.28   305.31     3.18     0.00   0.00726     71.4
     74671.3    2500.0   281.54   306.05     3.41     0.00   0.00653     70.2
     72442.7    2750.0   279.81   306.82     3.62     0.00   0.00587     69.0
     70266.6    3000.0   278.09   307.60     3.81     0.00   0.00527     67.8
     68142.3    3250.0   276.37   308.40     3.97     0.00   0.00473     66.6
     66069.0    3500.0   274.67   309.22     4.12     0.00   0.00423     65.3
     64045.8    3750.0   272.97   310.05     4.24     0.00   0.00378     64.0
     62072.1    4000.0   271.27   310.89     4.35     0.00   0.00337     62.7
     60147.0    4250.0   269.58   311.75     4.44     0.00   0.00300     61.4
     58269.6    4500.0   267.90   312.62     4.53     0.00   0.00267     60.1
     56439.2    4750.0   266.21   313.50     4.60     0.00   0.00237     58.8
     54654.9    5000.0   264.53   314.39     4.66     0.00   0.00210     57.4
     52916.0    5250.0   262.85   315.30     4.71     0.00   0.00185     56.1
     51221.5    5500.0   261.17   316.22     4.75     0.00   0.00163     54.7
     49570.7    5750.0   259.50   317.14     4.79     0.00   0.00143     53.3
     47962.7    6000.0   257.82   318.08     4.82     0.00   0.00126     51.9
     46396.7    6250.0   256.15   319.03     4.85     0.00   0.00110     50.5
     44872.0    6500.0   254.47   319.98     4.87     0.00   0.00096     49.1
     43387.6    6750.0   252.80   320.95     4.89     0.00   0.00083     47.6
     41942.9    7000.0   251.12   321.92     4.91     0.00   0.00072     46.2
     40537.0    7250.0   249.45   322.90     4.92     0.00   0.00062     44.7
     39169.2    7500.0   247.77   323.90     4.93     0.00   0.00054     43.2
     37838.6    7750.0   246.09   324.90     4.94     0.00   0.00046     41.8
     36544.5    8000.0   244.41   325.90     4.95     0.00   0.00039     40.3
     35286.2    8250.0   242.73   326.92     4.96     0.00   0.00033     38.8
     34062.8    8500.0   241.04   327.94     4.97     0.00   0.00028     37.3
     32873.7    8750.0   239.36   328.97     4.97     0.00   0.00024     35.8
     31718.1    9000.0   237.67   330.01     4.98     0.00   0.00020     34.3
     30595.3    9250.0   235.98   331.06     4.98     0.00   0.00017     32.8
     29504.5    9500.0   234.29   332.11     4.98     0.00   0.00014     31.3
     28445.1    9750.0   232.59   333.17     4.99     0.00   0.00012     29.8
     27416.4   10000.0   230.89   334.24     4.99     0.00   0.00010     28.3
     26417.7   10250.0   229.19   335.31     4.99     0.00   0.00008     26.8
     25448.3   10500.0   227.48   336.39     4.99     0.00   0.00006     25.3
     24507.6   10750.0   225.77   337.48     4.99     0.00   0.00005     23.8
     23594.9   11000.0   224.06   338.57     4.99     0.00   0.00004     22.3
     22709.5   11250.0   222.34   339.67     4.99     0.00   0.00003     20.8
     21850.9   11500.0   220.62   340.77     5.00     0.00   0.00003     19.4
     21018.4   11750.0   218.90   341.88     5.00     0.00   0.00002     17.9
     20211.4   12000.0   217.17   343.00     5.00     0.00   0.00002     16.5
     19432.4   12250.0   217.22   346.95     5.00     0.00   0.00001      9.9
     18683.6   12500.0   217.26   350.95     5.00     0.00   0.00001      9.5
     17963.8   12750.0   217.31   354.99     5.00     0.00   0.00001      9.1
     17272.0   13000.0   217.36   359.08     5.00     0.00   0.00001      8.7
     16606.9   13249.9   217.41   363.22     5.00     0.00   0.00001      8.3
     15967.5   13499.9   217.46   367.41     5.00     0.00   0.00001      7.9
     15353.0   13749.9   217.52   371.64     5.00     0.00   0.00001      7.6
     14762.2   13999.9   217.57   375.92     5.00     0.00   0.00001      7.2
     14194.3   14249.9   217.62   380.25     5.00     0.00   0.00001      6.9
     13648.3   14499.9   217.67   384.64     5.00     0.00   0.00001      6.6
     13123.5   14749.9   217.73   389.07     5.00     0.00   0.00001      6.3
     12619.0   14999.9   217.78   393.55     5.00     0.00   0.00001      6.0
     12134.0   15249.9   217.84   398.08     5.00     0.00   0.00001      5.7
     11667.8   15499.9   217.89   402.67     5.00     0.00   0.00001      5.5
     11219.6   15749.9   217.95   407.31     5.00     0.00   0.00001      5.2
     10788.7   15999.9   218.01   412.00     5.00     0.00   0.00001      5.0
     10374.5   16249.9   218.06   416.75     5.00     0.00   0.00001      4.8
      9976.2   16499.9   218.12   421.55     5.00     0.00   0.00001      4.6
      9593.4   16749.9   218.18   426.41     5.00     0.00   0.00001      4.4
      9225.4   16999.8   218.24   431.32     5.00     0.00   0.00001      4.2
      8871.5   17249.8   218.30   436.29     5.00     0.00   0.00001      4.0
      8531.4   17499.8   218.36   441.32     5.00     0.00   0.00001      3.8
      8204.4   17749.8   218.43   446.41     5.00     0.00   0.00001      3.6
      7889.9   17999.8   218.49   451.55     5.00     0.00   0.00001      3.5
      7587.7   18249.8   218.55   456.75     5.00     0.00   0.00001      3.3
      7297.1   18499.8   218.62   462.02     5.00     0.00   0.00001      3.1
      7017.7   18749.8   218.68   467.34     5.00     0.00   0.00001      3.0
      6749.0   18999.8   218.75   472.72     5.00     0.00   0.00001      2.9
      6490.8   19249.8   218.81   478.17     5.00     0.00   0.00001      2.7
      6242.5   19499.8   218.88   483.68     5.00     0.00   0.00001      2.6
      6003.7   19749.8   218.95   489.25     5.00     0.00   0.00001      2.5
      5774.2   19999.8   219.02   494.89     5.00     0.00   0.00001      2.4
      5553.5   20249.8   219.08   500.59     5.00     0.00   0.00001      2.3
      5341.3   20499.8   219.15   506.36     5.00     0.00   0.00001      2.2
      5137.3   20749.8   219.23   512.20     5.00     0.00   0.00001      2.1
      4941.1   20999.8   219.30   518.10     5.00     0.00   0.00001      2.0
      4752.5   21249.7   219.37   524.07     5.00     0.00   0.00001      1.9
      4571.2   21499.7   219.44   530.11     5.00     0.00   0.00001      1.8
      4396.8   21749.7   219.52   536.21     5.00     0.00   0.00001      1.7
      4229.1   21999.7   219.59   542.39     5.00     0.00   0.00001      1.6
      4067.9   22249.7   219.67   548.64     5.00     0.00   0.00001      1.5
      3912.9   22499.7   219.75   554.96     5.00     0.00   0.00001      1.5
      3763.8   22749.7   219.82   561.36     5.00     0.00   0.00001      1.4
      3620.5   22999.7   219.90   567.83     5.00     0.00   0.00001      1.3
      3482.7   23249.7   219.98   574.37     5.00     0.00   0.00001      1.3
      3350.1   23499.7   220.06   580.99     5.00     0.00   0.00001      1.2
      3222.7   23749.7   220.14   587.68     5.00     0.00   0.00001      1.2
      3100.1   23999.7   220.23   594.45     5.00     0.00   0.00001      1.1
      2982.3   24249.7   220.31   601.30     5.00     0.00   0.00001      1.0
      2869.0   24499.7   220.39   608.23     5.00     0.00   0.00001      1.0
      2760.0   24749.7   220.48   615.24     5.00     0.00   0.00001      1.0
      2655.2   24999.7   220.57   622.33     5.00     0.00   0.00001      0.9
      2554.4   25249.7   220.65   629.50     5.00     0.00   0.00001      0.9
      2457.5   25499.7   220.74   636.75     5.00     0.00   0.00001      0.8
      2364.3   25749.6   220.83   644.09     5.00     0.00   0.00001      0.8
      2274.7   25999.6   220.92   651.51     5.00     0.00   0.00001      0.7
      2188.5   26249.6   221.01   659.02     5.00     0.00   0.00001      0.7
      2105.6   26499.6   221.10   666.61     5.00     0.00   0.00001      0.7
      2025.9   26749.6   221.20   674.29     5.00     0.00   0.00001      0.6
      1949.2   26999.6   221.29   682.06     5.00     0.00   0.00001      0.6
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
Fsize = 20;
Wfilter = 20;

Fig = figure;
set(gca, 'FontSize', Fsize);
skew_sounding(Fig, P, T, TD, RH, WSPEED, WDIR, Wfilter);

%xlim([ -10 30 ]);
%ylim([ 500 1050 ]);

end