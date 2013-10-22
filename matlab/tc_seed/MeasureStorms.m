function [ ] = MeasureStorms(ConfigFile, InType, OutFile)
% MeasureStorms generate measurements of storm size, intensity
%

UseVt = strcmp(InType, 'vt');

[ Config ] = ReadConfig(ConfigFile);

UndefVal = Config.UndefVal;
Tdir = Config.TsavgDir;

if (strcmp(OutFile, 'stdout'))
  FileId = 1; % stdout
else
  FileId = fopen(OutFile, 'w');
end

% time interval for the "SS" phase
Tstart = 120;
Tend   = 140;

% compare only clean (C0100) and polluted (C2000) cases
Ccase = 'TCS_SD_C0100';
Pcase = 'TCS_SD_C2000';

KEvar = 'horiz_ke';
Vvar = 'max_azwind';
RMWvar = 'rmw';
R34KTvar = 'r34kt';
R50KTvar = 'r50kt';
R64KTvar = 'r64kt';

% Set file names based on tangential winds, or 10 m winds
if (UseVt)
  C_KEfile = sprintf('%s/%s_%s.h5', Tdir, KEvar, Ccase);
  P_KEfile = sprintf('%s/%s_%s.h5', Tdir, KEvar, Pcase);

  C_Vfile = sprintf('%s/%s_%s.h5', Tdir, Vvar, Ccase);
  P_Vfile = sprintf('%s/%s_%s.h5', Tdir, Vvar, Pcase);
else
  C_KEfile = sprintf('%s/%s_sp10_%s.h5', Tdir, KEvar, Ccase);
  P_KEfile = sprintf('%s/%s_sp10_%s.h5', Tdir, KEvar, Pcase);

  C_Vfile = sprintf('%s/%s_sp10_%s.h5', Tdir, Vvar, Ccase);
  P_Vfile = sprintf('%s/%s_sp10_%s.h5', Tdir, Vvar, Pcase);
end

fprintf(FileId,'Calculating percent change in storm stats\n');
if (UseVt)
  fprintf(FileId,'  Using TANGENTIAL Wind\n');
else
  fprintf(FileId,'  Using 10 METER Wind\n');
end
fprintf(FileId,'\n');

fprintf(FileId,'  Reading: %s, Dataset: %s\n', C_KEfile, KEvar);
fprintf(FileId,'  Reading: %s, Dataset: %s\n', P_KEfile, KEvar);
fprintf(FileId,'  Reading: %s, Dataset: %s\n', C_Vfile, Vvar);
fprintf(FileId,'  Reading: %s, Dataset: %s\n', P_Vfile, Vvar);
fprintf(FileId,'  Reading: %s, Dataset: %s\n', C_Vfile, RMWvar);
fprintf(FileId,'  Reading: %s, Dataset: %s\n', P_Vfile, RMWvar);
fprintf(FileId,'  Reading: %s, Dataset: %s\n', C_Vfile, R34KTvar);
fprintf(FileId,'  Reading: %s, Dataset: %s\n', P_Vfile, R34KTvar);
fprintf(FileId,'  Reading: %s, Dataset: %s\n', C_Vfile, R50KTvar);
fprintf(FileId,'  Reading: %s, Dataset: %s\n', P_Vfile, R50KTvar);
fprintf(FileId,'  Reading: %s, Dataset: %s\n', C_Vfile, R64KTvar);
fprintf(FileId,'  Reading: %s, Dataset: %s\n', P_Vfile, R64KTvar);
fprintf(FileId,'\n');

C_KE = squeeze(hdf5read(C_KEfile, KEvar));
P_KE = squeeze(hdf5read(P_KEfile, KEvar));
C_V = squeeze(hdf5read(C_Vfile, Vvar));
C_RMW = squeeze(hdf5read(C_Vfile, RMWvar));
C_R34KT = squeeze(hdf5read(C_Vfile, R34KTvar));
C_R50KT = squeeze(hdf5read(C_Vfile, R50KTvar));
C_R64KT = squeeze(hdf5read(C_Vfile, R64KTvar));
P_V = squeeze(hdf5read(P_Vfile, Vvar));
P_RMW = squeeze(hdf5read(P_Vfile, RMWvar));
P_R34KT = squeeze(hdf5read(P_Vfile, R34KTvar));
P_R50KT = squeeze(hdf5read(P_Vfile, R50KTvar));
P_R64KT = squeeze(hdf5read(P_Vfile, R64KTvar));

% Find time range for SS phase
T = squeeze(hdf5read(C_KEfile, 't_coords')) / 3600; % hr
T1 = find(T >= Tstart, 1, 'first');
T2 = find(T <= Tend,   1, 'last');

% calculate the percent changes
% Kinetic Energy
[ C_KE_LAST P_KE_LAST C_KE_AVG P_KE_AVG KE_LAST_CHG KE_AVG_CHG ]  = CalcPercentChange(C_KE, P_KE, T1, T2);

% Velocity
[ C_V_LAST P_V_LAST C_V_AVG P_V_AVG V_LAST_CHG V_AVG_CHG ]  = CalcPercentChange(C_V, P_V, T1, T2);

% radius of max winds, 34 kt, 50 kt and 64 kt
[ C_RMW_LAST P_RMW_LAST C_RMW_AVG P_RMW_AVG RMW_LAST_CHG RMW_AVG_CHG ]  = CalcPercentChange(C_RMW, P_RMW, T1, T2);
[ C_R34KT_LAST P_R34KT_LAST C_R34KT_AVG P_R34KT_AVG R34KT_LAST_CHG R34KT_AVG_CHG ]  = CalcPercentChange(C_R34KT, P_R34KT, T1, T2);
[ C_R50KT_LAST P_R50KT_LAST C_R50KT_AVG P_R50KT_AVG R50KT_LAST_CHG R50KT_AVG_CHG ]  = CalcPercentChange(C_R50KT, P_R50KT, T1, T2);
[ C_R64KT_LAST P_R64KT_LAST C_R64KT_AVG P_R64KT_AVG R64KT_LAST_CHG R64KT_AVG_CHG ]  = CalcPercentChange(C_R64KT, P_R64KT, T1, T2);

% print report
fprintf(FileId,'%15s %15s %15s %15s\n', 'Statistic', 'Clean', 'Polluted', '% Change');
fprintf(FileId,'\n');
fprintf(FileId,'%15s %15.2e %15.2e %15.2f\n', 'KE (last)', C_KE_LAST, P_KE_LAST, KE_LAST_CHG);
fprintf(FileId,'%15s %15.2e %15.2e %15.2f\n', 'KE (avg)', C_KE_AVG, P_KE_AVG, KE_AVG_CHG);
fprintf(FileId,'\n');
fprintf(FileId,'%15s %15.2e %15.2e %15.2f\n', 'V (last)', C_V_LAST, P_V_LAST, V_LAST_CHG);
fprintf(FileId,'%15s %15.2e %15.2e %15.2f\n', 'V (avg)', C_V_AVG, P_V_AVG, V_AVG_CHG);
fprintf(FileId,'\n');
fprintf(FileId,'%15s %15.2e %15.2e %15.2f\n', 'RMW (last)', C_RMW_LAST, P_RMW_LAST, RMW_LAST_CHG);
fprintf(FileId,'%15s %15.2e %15.2e %15.2f\n', 'RMW (avg)', C_RMW_AVG, P_RMW_AVG, RMW_AVG_CHG);
fprintf(FileId,'\n');
fprintf(FileId,'%15s %15.2e %15.2e %15.2f\n', 'R34KT (last)', C_R34KT_LAST, P_R34KT_LAST, R34KT_LAST_CHG);
fprintf(FileId,'%15s %15.2e %15.2e %15.2f\n', 'R34KT (avg)', C_R34KT_AVG, P_R34KT_AVG, R34KT_AVG_CHG);
fprintf(FileId,'\n');
fprintf(FileId,'%15s %15.2e %15.2e %15.2f\n', 'R50KT (last)', C_R50KT_LAST, P_R50KT_LAST, R50KT_LAST_CHG);
fprintf(FileId,'%15s %15.2e %15.2e %15.2f\n', 'R50KT (avg)', C_R50KT_AVG, P_R50KT_AVG, R50KT_AVG_CHG);
fprintf(FileId,'\n');
fprintf(FileId,'%15s %15.2e %15.2e %15.2f\n', 'R64KT (last)', C_R64KT_LAST, P_R64KT_LAST, R64KT_LAST_CHG);
fprintf(FileId,'%15s %15.2e %15.2e %15.2f\n', 'R64KT (avg)', C_R64KT_AVG, P_R64KT_AVG, R64KT_AVG_CHG);


% Close the file if not stdout nor stderr
if (FileId > 2)
  fclose(FileId);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CalcPercentChange()
%
function [ C_LAST P_LAST C_AVG P_AVG LAST_CHG AVG_CHG ] = CalcPercentChange(C, P, T1, T2)
% CalcPercentChange extract the last and average data values and calculate their percent change

C_LAST = C(end);
P_LAST = P(end);
C_AVG = mean(C(T1:T2));
P_AVG = mean(P(T1:T2));
LAST_CHG = ((P_LAST / C_LAST) - 1) * 100;
AVG_CHG  = ((P_AVG / C_AVG) - 1) * 100;

end
