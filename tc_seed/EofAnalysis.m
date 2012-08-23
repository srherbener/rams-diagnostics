% ATS655 Final project
%   EOF analysis on my own data
%
%  TC data
%  Time series of azimuthally averaged theta-e
%     Read in azavg hdf5 file --> (x,y,z,t)
%     Convert to 2D array for EOF routine --> (t,A) where A is linear storage of (x,z,t)
%     2D data: X-Z planes (50 x values, 18 z values --> 900 points per time step)
%     Data is organized: 2D array (145 rows, 900 columns)
%       rows --> time
%       columns --> 2D field
%          first 50 entries are x = 1:50, z = 1
%          next 50 entries are x = 1:50, z = 2
%          etc.
%
%  Time series of storm intensity metric
%     1D data: (145 rows, 1 column)
%       rows --> time
%       columns --> metric value
%
%  Data sets for two experiments varying CCN concentrations: 100/cc and
%  8000/cc

clear;

% Read in the HDF5 data --> (x,y,z,t)
H5_CCN_CNTL  = hdf5read('AzAveragedData/ccn_conc_TCS_CNTL.h5','ccn_conc');
H5_CCN_C0100 = hdf5read('AzAveragedData/ccn_conc_TCS_GN_C0100.h5','ccn_conc');
H5_CCN_C0500 = hdf5read('AzAveragedData/ccn_conc_TCS_GN_C0500.h5','ccn_conc');
H5_CCN_C1000 = hdf5read('AzAveragedData/ccn_conc_TCS_GN_C1000.h5','ccn_conc');
H5_CCN_C2000 = hdf5read('AzAveragedData/ccn_conc_TCS_GN_C2000.h5','ccn_conc');

H5_TE_CNTL  = hdf5read('AzAveragedData/theta_e_TCS_CNTL.h5','theta_e');
H5_TE_C0100 = hdf5read('AzAveragedData/theta_e_TCS_GN_C0100.h5','theta_e');
H5_TE_C0500 = hdf5read('AzAveragedData/theta_e_TCS_GN_C0500.h5','theta_e');
H5_TE_C1000 = hdf5read('AzAveragedData/theta_e_TCS_GN_C1000.h5','theta_e');
H5_TE_C2000 = hdf5read('AzAveragedData/theta_e_TCS_GN_C2000.h5','theta_e');

H5_XCOORDS = hdf5read('AzAveragedData/theta_e_TCS_CNTL.h5', 'x_coords');
H5_ZCOORDS = hdf5read('AzAveragedData/theta_e_TCS_CNTL.h5', 'z_coords');

% Figure out the selection of the (x,y,z,t) data
[ Nx, Ny, Nz, Nt ] = size(H5_TE_CNTL);
X1 = 4; % exclude the 4km, 8km, 12km radius points --> exclude warm core
X2 = Nx;
Y1 = 1;
Y2 = Ny;
Z1 = find(H5_ZCOORDS > 0, 1);           % look at low level theta_e
Z2 = find(H5_ZCOORDS < 4000, 1, 'last');
CCN_Z2 = find(H5_ZCOORDS < 10000, 1, 'last');
T1 = 1;
T2 = Nt;

% Coordinate values for the plots
Radii = H5_XCOORDS(X1:X2) / 1000; % convert from m to km
Heights = H5_ZCOORDS(Z1:Z2);
CCN_Heights = H5_ZCOORDS(Z1:CCN_Z2);

% Convert to the 2D format that EOF analysis wants --> (t, A)
TE_CNTL  = Xyzt2EofArray(H5_TE_CNTL, X1, X2, Y1, Y2, Z1, Z2, T1, T2);
TE_C0100 = Xyzt2EofArray(H5_TE_C0100, X1, X2, Y1, Y2, Z1, Z2, T1, T2);
TE_C0500 = Xyzt2EofArray(H5_TE_C0500, X1, X2, Y1, Y2, Z1, Z2, T1, T2);
TE_C1000 = Xyzt2EofArray(H5_TE_C1000, X1, X2, Y1, Y2, Z1, Z2, T1, T2);
TE_C2000 = Xyzt2EofArray(H5_TE_C2000, X1, X2, Y1, Y2, Z1, Z2, T1, T2);

CCN_CNTL  = Xyzt2EofArray(H5_CCN_CNTL, X1, X2, Y1, Y2, Z1, CCN_Z2, T1, T2);
CCN_C0100 = Xyzt2EofArray(H5_CCN_C0100, X1, X2, Y1, Y2, Z1, CCN_Z2, T1, T2);
CCN_C0500 = Xyzt2EofArray(H5_CCN_C0500, X1, X2, Y1, Y2, Z1, CCN_Z2, T1, T2);
CCN_C1000 = Xyzt2EofArray(H5_CCN_C1000, X1, X2, Y1, Y2, Z1, CCN_Z2, T1, T2);
CCN_C2000 = Xyzt2EofArray(H5_CCN_C2000, X1, X2, Y1, Y2, Z1, CCN_Z2, T1, T2);

% convert the undefined value (-999) to nan
TE_CNTL(TE_CNTL == -999) = nan;
TE_C0100(TE_C0100 == -999) = nan;
TE_C0500(TE_C0500 == -999) = nan;
TE_C1000(TE_C1000 == -999) = nan;
TE_C2000(TE_C2000 == -999) = nan;

CCN_CNTL(CCN_CNTL == -999) = nan;
CCN_C0100(CCN_C0100 == -999) = nan;
CCN_C0500(CCN_C0500 == -999) = nan;
CCN_C1000(CCN_C1000 == -999) = nan;
CCN_C2000(CCN_C2000 == -999) = nan;

% Form the difference between the experiments and the control, and run the
% anlysis on the variance of the differences.
ThetaE_C0100 = TE_C0100 - TE_CNTL;
ThetaE_C0500 = TE_C0500 - TE_CNTL;
ThetaE_C1000 = TE_C1000 - TE_CNTL;
ThetaE_C2000 = TE_C2000 - TE_CNTL;

Ccn_C0100 = CCN_C0100 - CCN_CNTL;
Ccn_C0500 = CCN_C0500 - CCN_CNTL;
Ccn_C1000 = CCN_C1000 - CCN_CNTL;
Ccn_C2000 = CCN_C2000 - CCN_CNTL;

% This might be bogus
Ccn_C0100(isnan(Ccn_C0100)) = 0;
Ccn_C0500(isnan(Ccn_C0500)) = 0;
Ccn_C1000(isnan(Ccn_C1000)) = 0;
Ccn_C2000(isnan(Ccn_C2000)) = 0;

% Run EOF on the temperature and theta_e data
[ EOF_TE_C0100, PC_TE_C0100, EV_TE_C0100 ] = EofNan(ThetaE_C0100,1);
[ EOF_TE_C0500, PC_TE_C0500, EV_TE_C0500 ] = EofNan(ThetaE_C0500,1);
[ EOF_TE_C1000, PC_TE_C1000, EV_TE_C1000 ] = EofNan(ThetaE_C1000,1);
[ EOF_TE_C2000, PC_TE_C2000, EV_TE_C2000 ] = EofNan(ThetaE_C2000,1);

[ EOF_CCN_C0100, PC_CCN_C0100, EV_CCN_C0100 ] = EofNan(Ccn_C0100,1);
[ EOF_CCN_C0500, PC_CCN_C0500, EV_CCN_C0500 ] = EofNan(Ccn_C0500,1);
[ EOF_CCN_C1000, PC_CCN_C1000, EV_CCN_C1000 ] = EofNan(Ccn_C1000,1);
[ EOF_CCN_C2000, PC_CCN_C2000, EV_CCN_C2000 ] = EofNan(Ccn_C2000,1);

%%% % Sweep lag correlations between PCs of temp and theta-e versus the storm
%%% % intensity metric.
%%% Lags = (1:50);
%%% [ LagCors_T, LagRcoeffs_T, LagRyints_T ] = SweepLagCor(StmInt,PC_T(:,1),Lags);
%%% [ LagCors_ThetaE, LagRcoeffs_ThetaE, LagRyints_ThetaE ] = SweepLagCor(StmInt,PC_ThetaE(:,1),Lags);
%%% [ LagCors_CCN, LagRcoeffs_CCN, LagRyints_CCN ] = SweepLagCor(StmInt,PC_CCN(:,1),Lags);
%%% 
%%% % Plots
%%% FigSi =figure;
%%% Tsteps = (1:length(SI_C100));
%%% plot(Tsteps, SI_C100);
%%% line(Tsteps, SI_C8000, 'Color', 'r');
%%% line(Tsteps, StmInt, 'Color', 'k', 'Linestyle', '--');
%%% title('Storm Intensity');
%%% xlabel('Timestep (each interval is 10 minutes)');
%%% ylabel('Intensity Metric');
%%% legend('CCN concentration: 100/cc', 'CCN concentration: 8000/cc', 'Difference (100/cc case - 8000/cc case)', 'Location', 'NorthWest');
%%% saveas(FigSi, 'SI_C100_C800.jpg');
%%% 

% THETA-E
EofNum = 1;

Clim = 0.05;
Cinc = Clim / 20;
Clevs = (-Clim:Cinc:Clim);
Cbounds = [ -Clim Clim ];

CCN_Clim = 0.1;
CCN_Cinc = CCN_Clim / 40;
CCN_Clevs = (-CCN_Clim:CCN_Cinc:CCN_Clim);
CCN_Cbounds = [ -CCN_Clim CCN_Clim ];

% Look at lower levels in the Theta-E EOF data
% Needed to include 0 - 10000m to get statistical significance of the first EOF, but all
% interesting activity with Theta-E is near the surface.
EZ = find(Heights < 1200, 1, 'last');
EP = length(Radii) * EZ;

Plot2dEofPc(EOF_TE_C0100(1:EP,:), PC_TE_C0100, EofNum, Radii, Heights(1:EZ), Clevs, Cbounds, 'Theta-e Difference', 'C0100', 'Radius (km)', 'Height (m)', 'Theta-e (K)', 'plots/EOF_ThetaE_C0100.jpg')
Plot2dEofPc(EOF_TE_C0500(1:EP,:), PC_TE_C0500, EofNum, Radii, Heights(1:EZ), Clevs, Cbounds, 'Theta-e Difference', 'C0500', 'Radius (km)', 'Height (m)', 'Theta-e (K)', 'plots/EOF_ThetaE_C0500.jpg')
Plot2dEofPc(EOF_TE_C1000(1:EP,:), PC_TE_C1000, EofNum, Radii, Heights(1:EZ), Clevs, Cbounds, 'Theta-e Difference', 'C1000', 'Radius (km)', 'Height (m)', 'Theta-e (K)', 'plots/EOF_ThetaE_C1000.jpg')
Plot2dEofPc(EOF_TE_C2000(1:EP,:), PC_TE_C2000, EofNum, Radii, Heights(1:EZ), Clevs, Cbounds, 'Theta-e Difference', 'C2000', 'Radius (km)', 'Height (m)', 'Theta-e (K)', 'plots/EOF_ThetaE_C2000.jpg')

% CCN
Plot2dEofPc(EOF_CCN_C0100, PC_CCN_C0100, EofNum, Radii, CCN_Heights, CCN_Clevs, CCN_Cbounds, 'CCN Concentration Difference', 'C0100', 'Radius (km)', 'Height (m)', 'CCN Concentration (#/cc)', 'plots/EOF_CCN_C0100.jpg')
Plot2dEofPc(EOF_CCN_C0500, PC_CCN_C0500, EofNum, Radii, CCN_Heights, CCN_Clevs, CCN_Cbounds, 'CCN Concentration Difference', 'C0500', 'Radius (km)', 'Height (m)', 'CCN Concentration (#/cc)', 'plots/EOF_CCN_C0500.jpg')
Plot2dEofPc(EOF_CCN_C1000, PC_CCN_C1000, EofNum, Radii, CCN_Heights, CCN_Clevs, CCN_Cbounds, 'CCN Concentration Difference', 'C1000', 'Radius (km)', 'Height (m)', 'CCN Concentration (#/cc)', 'plots/EOF_CCN_C1000.jpg')
Plot2dEofPc(EOF_CCN_C2000, PC_CCN_C2000, EofNum, Radii, CCN_Heights, CCN_Clevs, CCN_Cbounds, 'CCN Concentration Difference', 'C2000', 'Radius (km)', 'Height (m)', 'CCN Concentration (#/cc)', 'plots/EOF_CCN_C2000.jpg')

% Generate and plot Eigenvalue spectra
% Use N* = 20 (timestep == 1hr, 121 time steps --> ~6hr to get rid of persistence)
Nstar = 20;
NumEv = 10;

GenPlotEigSpectrum(EV_TE_C0100, Nstar, NumEv, 'Theta-e Difference', 'C0100', 'plots/ES_TE_C0100.jpg')
GenPlotEigSpectrum(EV_TE_C0500, Nstar, NumEv, 'Theta-e Difference', 'C0500', 'plots/ES_TE_C0500.jpg')
GenPlotEigSpectrum(EV_TE_C1000, Nstar, NumEv, 'Theta-e Difference', 'C1000', 'plots/ES_TE_C1000.jpg')
GenPlotEigSpectrum(EV_TE_C2000, Nstar, NumEv, 'Theta-e Difference', 'C2000', 'plots/ES_TE_C2000.jpg')

GenPlotEigSpectrum(EV_CCN_C0100, Nstar, NumEv, 'CCN Concentration Difference', 'C0100', 'plots/ES_CCN_C0100.jpg')
GenPlotEigSpectrum(EV_CCN_C0500, Nstar, NumEv, 'CCN Concentration Difference', 'C0500', 'plots/ES_CCN_C0500.jpg')
GenPlotEigSpectrum(EV_CCN_C1000, Nstar, NumEv, 'CCN Concentration Difference', 'C1000', 'plots/ES_CCN_C1000.jpg')
GenPlotEigSpectrum(EV_CCN_C2000, Nstar, NumEv, 'CCN Concentration Difference', 'C2000', 'plots/ES_CCN_C2000.jpg')

%%% FigLagCorT = figure;
%%% plot(Lags,LagCors_T);
%%% title('Lag Correlations: Temperature Difference vs Storm Intensity Difference');
%%% xlabel('Lag(each interval is 10 minutes)');
%%% ylabel('Correlation');
%%% saveas(FigLagCorT,'LagCor_T.jpg');
%%% 
%%% FigLagCorThetaE = figure;
%%% plot(Lags,LagCors_ThetaE);
%%% title('Lag Correlations: Theta-e Difference vs Storm Intensity Difference');
%%% xlabel('Lag(each interval is 10 minutes)');
%%% ylabel('Correlation');
%%% saveas(FigLagCorThetaE,'LagCor_ThetaE.jpg');
%%% 
%%% FigLagCorCcn = figure;
%%% plot(Lags,LagCors_CCN);
%%% title('Lag Correlations: CCN Concentration Difference vs Storm Intensity Difference');
%%% xlabel('Lag(each interval is 10 minutes)');
%%% ylabel('Correlation');
%%% saveas(FigLagCorCcn,'LagCor_CCN.jpg');
%%% 
