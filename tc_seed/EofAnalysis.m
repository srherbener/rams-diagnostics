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
H5_TE_CNTL  = hdf5read('AzAveragedData/theta_e_TCS_CNTL.h5','theta_e');
H5_TE_C0100 = hdf5read('AzAveragedData/theta_e_TCS_GN_C0100.h5','theta_e');
H5_TE_C0500 = hdf5read('AzAveragedData/theta_e_TCS_GN_C0500.h5','theta_e');
H5_TE_C1000 = hdf5read('AzAveragedData/theta_e_TCS_GN_C1000.h5','theta_e');
H5_TE_C2000 = hdf5read('AzAveragedData/theta_e_TCS_GN_C2000.h5','theta_e');

H5_XCOORDS = hdf5read('AzAveragedData/theta_e_TCS_CNTL.h5', 'x_coords');
H5_ZCOORDS = hdf5read('AzAveragedData/theta_e_TCS_CNTL.h5', 'z_coords');

% Figure out the selection of the (x,y,z,t) data
[ Nx, Ny, Nz, Nt ] = size(H5_TE_CNTL);
X1 = 1;
X2 = Nx;
Y1 = 1;
Y2 = Ny;
Z1 = find(H5_ZCOORDS > 0, 1);           % look at low level theta_e
Z2 = find(H5_ZCOORDS < 2000, 1, 'last');
T1 = 1;
T2 = Nt;

% Coordinate values for the plots
Radii = H5_XCOORDS(X1:X2);
Heights = H5_ZCOORDS(Z1:Z2);

% Convert to the 2D format that EOF analysis wants --> (t, A)
TE_CNTL  = Xyzt2EofArray(H5_TE_CNTL, X1, X2, Y1, Y2, Z1, Z2, T1, T2);
TE_C0100 = Xyzt2EofArray(H5_TE_C0100, X1, X2, Y1, Y2, Z1, Z2, T1, T2);
TE_C0500 = Xyzt2EofArray(H5_TE_C0500, X1, X2, Y1, Y2, Z1, Z2, T1, T2);
TE_C1000 = Xyzt2EofArray(H5_TE_C1000, X1, X2, Y1, Y2, Z1, Z2, T1, T2);
TE_C2000 = Xyzt2EofArray(H5_TE_C2000, X1, X2, Y1, Y2, Z1, Z2, T1, T2);

% Form the difference between the control and experiments and run the
% anlysis on the variance of the differences.
ThetaE_C0100 = TE_CNTL - TE_C0100;
ThetaE_C0500 = TE_CNTL - TE_C0500;
ThetaE_C1000 = TE_CNTL - TE_C1000;
ThetaE_C2000 = TE_CNTL - TE_C2000;

% Run EOF on the temperature and theta_e data
[ EOF_TE_C0100, PC_TE_C0100, EV_TE_C0100 ] = EofNan(ThetaE_C0100,1);
[ EOF_TE_C0500, PC_TE_C0500, EV_TE_C0500 ] = EofNan(ThetaE_C0500,1);
[ EOF_TE_C1000, PC_TE_C1000, EV_TE_C1000 ] = EofNan(ThetaE_C1000,1);
[ EOF_TE_C2000, PC_TE_C2000, EV_TE_C2000 ] = EofNan(ThetaE_C2000,1);


%%% % Form eigenvalue spectra
%%% % Use N* = 10 (guessing for now)
%%% Nstar = 10;
%%% %[ Lambdas_T_C100, VarExpl_T_C100, Err_T_C100 ] = EigenSpectrum(EV_T_C100, Nstar);
%%% %[ Lambdas_ThetaE_C100, VarExpl_ThetaE_C100, Err_ThetaE_C100 ] = EigenSpectrum(EV_ThetaE_C100, Nstar);
%%% %[ Lambdas_T_C8000, VarExpl_T_C8000, Err_T_C8000 ] = EigenSpectrum(EV_T_C8000, Nstar);
%%% %[ Lambdas_ThetaE_C8000, VarExpl_ThetaE_C8000, Err_ThetaE_C8000 ] = EigenSpectrum(EV_ThetaE_C8000, Nstar);
%%% [ Lambda_T, VarExpl_T, Err_T ] = EigenSpectrum(EV_T, Nstar);
%%% [ Lambda_ThetaE, VarExpl_ThetaE, Err_ThetaE ] = EigenSpectrum(EV_ThetaE, Nstar);
%%% [ Lambda_CCN, VarExpl_CCN, Err_CCN ] = EigenSpectrum(EV_CCN, Nstar);
%%% 
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
%%% FigEofT = figure;
%%% subplot(2,1,1);
%%% Plot2dMap(FigEofT,Radii,Heights,EOF_T(:,1)','Radius (km)', 'Height (m)', 'Temperature Difference EOF1');
%%% subplot(2,1,2);
%%% plot(PC_T(:,1)');
%%% title('Temperature Difference PC1');
%%% xlabel('Timestep');
%%% ylabel('Temperature (\circC)');
%%% saveas(FigEofT, 'EOF_Temp.jpg');
%%% 

FigEofTE_C0100 = figure;
subplot(2,1,1);
Plot2dMap(FigEofTE_C0100,Radii,Heights,EOF_TE_C0100(:,1)','Radius (km)', 'Height (m)', 'Theta-e Difference EOF1: C0100');
subplot(2,1,2);
plot(PC_TE_C0100(:,1)');
title('Theta-e Difference PC1');
xlabel('Timestep');
ylabel('Theta-e (K)');
saveas(FigEofThetaE, 'EOF_ThetaE_C0100.jpg');

%%% 
%%% FigEofCcn = figure;
%%% subplot(2,1,1);
%%% Plot2dMap(FigEofCcn,Radii,Heights,EOF_CCN(:,1)','Radius (km)', 'Height (m)', 'CCN Concentration Difference EOF1');
%%% subplot(2,1,2);
%%% plot(PC_CCN(:,1)');
%%% title('CCN Concentration Difference PC1');
%%% xlabel('Timestep');
%%% ylabel('CCN Concentration (#/cc)');
%%% saveas(FigEofCcn, 'EOF_CCN.jpg');
%%% 
%%% FigEsT = figure;
%%% errorbar(VarExpl_T(1:10), Err_T(1:10)); % just show the first 10
%%% title('First 10 Eigenvalues of Temperature Difference covariance (% Variance Explained)');
%%% xlabel('Eigenvalue number');
%%% ylabel('% Variance Explained');
%%% saveas(FigEsT, 'ES_Temp.jpg');
%%% 
%%% FigEsThetaE = figure;
%%% errorbar(VarExpl_ThetaE(1:10), Err_ThetaE(1:10)); % just show the first 10
%%% title('First 10 Eigenvalues of Theta-e Difference covariance (% Variance Explained)');
%%% xlabel('Eigenvalue number');
%%% ylabel('% Variance Explained');
%%% saveas(FigEsThetaE, 'ES_ThetaE.jpg');
%%% 
%%% FigEsCcn = figure;
%%% errorbar(VarExpl_CCN(1:10), Err_CCN(1:10)); % just show the first 10
%%% title('First 10 Eigenvalues of CCN Concentration Difference covariance (% Variance Explained)');
%%% xlabel('Eigenvalue number');
%%% ylabel('% Variance Explained');
%%% saveas(FigEsCcn, 'ES_CCN.jpg');
%%% 
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
