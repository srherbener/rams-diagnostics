% ATS655 Final project
%   EOF analysis on my own data
%
%  TC data
%  Time series of azimuthally averaged theta-e, temp (degrees C)
%     2D data: X-Z planes (50 x values, 18 z values --> 900 points per time
%     step)
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
path(path,'../Matlab');

% Read in the data
Raw_T_C100 = importdata('TC_SEED_C0100_tempc.txt');
Raw_ThetaE_C100 = importdata('TC_SEED_C0100_theta_e.txt');
Raw_CCN_C100 = importdata('TC_SEED_C0100_ccn_conc.txt');
SI_C100 = importdata('TC_SEED_C0100_ts_storm_int.txt');

Raw_T_C8000 = importdata('TC_SEED_C8000_tempc.txt');
Raw_ThetaE_C8000 = importdata('TC_SEED_C8000_theta_e.txt');
Raw_CCN_C8000 = importdata('TC_SEED_C8000_ccn_conc.txt');
SI_C8000 = importdata('TC_SEED_C8000_ts_storm_int.txt');

Radii = importdata('TC_SEED_Xvals.txt');
Heights = importdata('TC_SEED_Yvals.txt');

% Remove the time (column) means from the temperature and theta_e data.
%T_C100 = RemoveMeanNan(Raw_T_C100,1);
%ThetaE_C100 = RemoveMeanNan(Raw_ThetaE_C100,1);
%T_C8000 = RemoveMeanNan(Raw_T_C8000,1);
%ThetaE_C8000 = RemoveMeanNan(Raw_ThetaE_C8000,1);
%T_C100 = Raw_T_C100;
%ThetaE_C100 = Raw_ThetaE_C100;
%T_C8000 = Raw_T_C8000;
%ThetaE_C8000 = Raw_ThetaE_C8000;

% Form the difference between the C100 and C8000 experiments and run the
% anlysis on the variance of the differences.
Temp = Raw_T_C100 - Raw_T_C8000;
ThetaE = Raw_ThetaE_C100 - Raw_ThetaE_C8000;
CCN = Raw_CCN_C100 - Raw_CCN_C8000;
StmInt = SI_C100 - SI_C8000;

% Run EOF on the temperature and theta_e data
%[ EOF_T_C100, PC_T_C100, EV_T_C100 ] = EofNan(T_C100, 1);
%[ EOF_ThetaE_C100, PC_ThetaE_C100, EV_ThetaE_C100 ] = EofNan(ThetaE_C100, 1);
%[ EOF_T_C8000, PC_T_C8000, EV_T_C8000 ] = EofNan(T_C8000, 1);
%[ EOF_ThetaE_C8000, PC_ThetaE_C8000, EV_ThetaE_C8000 ] = EofNan(ThetaE_C8000, 1);
[ EOF_T, PC_T, EV_T ] = EofNan(Temp,1);
[ EOF_ThetaE, PC_ThetaE, EV_ThetaE ] = EofNan(ThetaE,1);
[ EOF_CCN, PC_CCN, EV_CCN ] = EofNan(CCN,1);

% Form eigenvalue spectra
% Use N* = 10 (guessing for now)
Nstar = 10;
%[ Lambdas_T_C100, VarExpl_T_C100, Err_T_C100 ] = EigenSpectrum(EV_T_C100, Nstar);
%[ Lambdas_ThetaE_C100, VarExpl_ThetaE_C100, Err_ThetaE_C100 ] = EigenSpectrum(EV_ThetaE_C100, Nstar);
%[ Lambdas_T_C8000, VarExpl_T_C8000, Err_T_C8000 ] = EigenSpectrum(EV_T_C8000, Nstar);
%[ Lambdas_ThetaE_C8000, VarExpl_ThetaE_C8000, Err_ThetaE_C8000 ] = EigenSpectrum(EV_ThetaE_C8000, Nstar);
[ Lambda_T, VarExpl_T, Err_T ] = EigenSpectrum(EV_T, Nstar);
[ Lambda_ThetaE, VarExpl_ThetaE, Err_ThetaE ] = EigenSpectrum(EV_ThetaE, Nstar);
[ Lambda_CCN, VarExpl_CCN, Err_CCN ] = EigenSpectrum(EV_CCN, Nstar);

% Sweep lag correlations between PCs of temp and theta-e versus the storm
% intensity metric.
Lags = (1:50);
[ LagCors_T, LagRcoeffs_T, LagRyints_T ] = SweepLagCor(StmInt,PC_T(:,1),Lags);
[ LagCors_ThetaE, LagRcoeffs_ThetaE, LagRyints_ThetaE ] = SweepLagCor(StmInt,PC_ThetaE(:,1),Lags);
[ LagCors_CCN, LagRcoeffs_CCN, LagRyints_CCN ] = SweepLagCor(StmInt,PC_CCN(:,1),Lags);

% Plots
FigSi =figure;
Tsteps = (1:length(SI_C100));
plot(Tsteps, SI_C100);
line(Tsteps, SI_C8000, 'Color', 'r');
line(Tsteps, StmInt, 'Color', 'k', 'Linestyle', '--');
title('Storm Intensity');
xlabel('Timestep (each interval is 10 minutes)');
ylabel('Intensity Metric');
legend('CCN concentration: 100/cc', 'CCN concentration: 8000/cc', 'Difference (100/cc case - 8000/cc case)', 'Location', 'NorthWest');
saveas(FigSi, 'SI_C100_C800.jpg');

FigEofT = figure;
subplot(2,1,1);
Plot2dMap(FigEofT,Radii,Heights,EOF_T(:,1)','Radius (km)', 'Height (m)', 'Temperature Difference EOF1');
subplot(2,1,2);
plot(PC_T(:,1)');
title('Temperature Difference PC1');
xlabel('Timestep');
ylabel('Temperature (\circC)');
saveas(FigEofT, 'EOF_Temp.jpg');

FigEofThetaE = figure;
subplot(2,1,1);
Plot2dMap(FigEofThetaE,Radii,Heights,EOF_ThetaE(:,1)','Radius (km)', 'Height (m)', 'Theta-e Difference EOF1');
subplot(2,1,2);
plot(PC_ThetaE(:,1)');
title('Theta-e Difference PC1');
xlabel('Timestep');
ylabel('Theta-e (K)');
saveas(FigEofThetaE, 'EOF_ThetaE.jpg');

FigEofCcn = figure;
subplot(2,1,1);
Plot2dMap(FigEofCcn,Radii,Heights,EOF_CCN(:,1)','Radius (km)', 'Height (m)', 'CCN Concentration Difference EOF1');
subplot(2,1,2);
plot(PC_CCN(:,1)');
title('CCN Concentration Difference PC1');
xlabel('Timestep');
ylabel('CCN Concentration (#/cc)');
saveas(FigEofCcn, 'EOF_CCN.jpg');

FigEsT = figure;
errorbar(VarExpl_T(1:10), Err_T(1:10)); % just show the first 10
title('First 10 Eigenvalues of Temperature Difference covariance (% Variance Explained)');
xlabel('Eigenvalue number');
ylabel('% Variance Explained');
saveas(FigEsT, 'ES_Temp.jpg');

FigEsThetaE = figure;
errorbar(VarExpl_ThetaE(1:10), Err_ThetaE(1:10)); % just show the first 10
title('First 10 Eigenvalues of Theta-e Difference covariance (% Variance Explained)');
xlabel('Eigenvalue number');
ylabel('% Variance Explained');
saveas(FigEsThetaE, 'ES_ThetaE.jpg');

FigEsCcn = figure;
errorbar(VarExpl_CCN(1:10), Err_CCN(1:10)); % just show the first 10
title('First 10 Eigenvalues of CCN Concentration Difference covariance (% Variance Explained)');
xlabel('Eigenvalue number');
ylabel('% Variance Explained');
saveas(FigEsCcn, 'ES_CCN.jpg');

FigLagCorT = figure;
plot(Lags,LagCors_T);
title('Lag Correlations: Temperature Difference vs Storm Intensity Difference');
xlabel('Lag(each interval is 10 minutes)');
ylabel('Correlation');
saveas(FigLagCorT,'LagCor_T.jpg');

FigLagCorThetaE = figure;
plot(Lags,LagCors_ThetaE);
title('Lag Correlations: Theta-e Difference vs Storm Intensity Difference');
xlabel('Lag(each interval is 10 minutes)');
ylabel('Correlation');
saveas(FigLagCorThetaE,'LagCor_ThetaE.jpg');

FigLagCorCcn = figure;
plot(Lags,LagCors_CCN);
title('Lag Correlations: CCN Concentration Difference vs Storm Intensity Difference');
xlabel('Lag(each interval is 10 minutes)');
ylabel('Correlation');
saveas(FigLagCorCcn,'LagCor_CCN.jpg');

