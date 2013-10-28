function [ ] = MeasureStorms(ConfigFile, InType, OutFile)
% MeasureStorms generate measurements of storm size, intensity
%

UseVt = strcmp(InType, 'vt');

[ Config ] = ReadConfig(ConfigFile);

UndefVal = Config.UndefVal;
Adir = Config.AzavgDir;
Tdir = Config.TsavgDir;
ControlCase = Config.ControlCase;

if (strcmp(OutFile, 'stdout'))
  FileId = 1; % stdout
else
  FileId = fopen(OutFile, 'w');
end

% time interval for the "SS" phase
Tstart = 120;
Tend   = 140;

KEvar = 'horiz_ke';
Vvar = 'max_azwind';
RMWvar = 'rmw';
R34KTvar = 'r34kt';
R50KTvar = 'r50kt';
R64KTvar = 'r64kt';
Pvar = 'sea_press';

fprintf(FileId,'Measuring storm structural quantities\n');
if (UseVt)
  fprintf(FileId,'  Using TANGENTIAL Wind\n');
else
  fprintf(FileId,'  Using 10 METER Wind\n');
end
fprintf(FileId,'\n');

for ic = 1:length(Config.Cases)
  Case = Config.Cases(ic).Cname;
  CASES{ic} = Case;

  % Set file names based on tangential winds, or 10 m winds
  if (UseVt)
    KEfile = sprintf('%s/%s_%s.h5', Tdir, KEvar, Case);
    Vfile = sprintf('%s/%s_%s.h5', Tdir, Vvar, Case);
  else
    KEfile = sprintf('%s/%s_sp10_%s.h5', Tdir, KEvar, Case);
    Vfile = sprintf('%s/%s_sp10_%s.h5', Tdir, Vvar, Case);
  end
  Pfile = sprintf('%s/%s_%s.h5', Adir, Pvar, Case);

  fprintf(FileId,'  Reading: %s, Dataset: %s\n', KEfile, KEvar);
  fprintf(FileId,'  Reading: %s, Dataset: %s\n', Vfile, Vvar);
  fprintf(FileId,'  Reading: %s, Dataset: %s\n', Vfile, RMWvar);
  fprintf(FileId,'  Reading: %s, Dataset: %s\n', Vfile, R64KTvar);
  fprintf(FileId,'  Reading: %s, Dataset: %s\n', Vfile, R50KTvar);
  fprintf(FileId,'  Reading: %s, Dataset: %s\n', Vfile, R34KTvar);
  fprintf(FileId,'  Reading: %s, Dataset: %s\n', Pfile, Pvar);
  fprintf(FileId,'\n');

  KE    = squeeze(hdf5read(KEfile, KEvar));
  V     = squeeze(hdf5read(Vfile, Vvar));
  RMW   = squeeze(hdf5read(Vfile, RMWvar));
  R64KT = squeeze(hdf5read(Vfile, R64KTvar));
  R50KT = squeeze(hdf5read(Vfile, R50KTvar));
  R34KT = squeeze(hdf5read(Vfile, R34KTvar));
  P     = squeeze(hdf5read(Pfile, Pvar));

  % Find time range for SS phase
  T = squeeze(hdf5read(KEfile, 't_coords')) / 3600; % hr
  T1 = find(T >= Tstart, 1, 'first');
  T2 = find(T <= Tend,   1, 'last');

  % record stats for the last time value, and for the average between T1 and T2
  KE_LAST(ic)    = KE(end);
  KE_AVG(ic)     = mean(KE(T1:T2));
  V_LAST(ic)     = V(end);
  V_AVG(ic)      = mean(V(T1:T2));
  RMW_LAST(ic)   = RMW(end);
  RMW_AVG(ic)    = mean(RMW(T1:T2));
  R64KT_LAST(ic) = R64KT(end);
  R64KT_AVG(ic)  = mean(R64KT(T1:T2));
  R50KT_LAST(ic) = R50KT(end);
  R50KT_AVG(ic)  = mean(R50KT(T1:T2));
  R34KT_LAST(ic) = R34KT(end);
  R34KT_AVG(ic)  = mean(R34KT(T1:T2));

  % MSP will be organized as (r,t)
  % first find the min value along r, then the last and mean time values
  MIN_P = squeeze(min(P,[],1));
  MSP_LAST(ic)   = MIN_P(end);
  MSP_AVG(ic)    = mean(MIN_P(T1:T2));
end

% print report
fprintf(FileId, 'Quantities taken from last time step:\n');
fprintf(FileId,'\n');
fprintf(FileId, '%15s %10s %10s %10s %10s %10s %10s %10s\n', 'Case', 'IKE', 'RMW', 'R64KT', 'R50KT', 'R34KT', 'V10', 'MSP');
fprintf(FileId,'\n');
for ic = 1:length(CASES)
  fprintf(FileId,'%15s %10.2e %10.2e %10.2e %10.2e %10.2e %10.2e %10.2e\n', ...
    CASES{ic}, KE_LAST(ic), RMW_LAST(ic), R64KT_LAST(ic), R50KT_LAST(ic), R34KT_LAST(ic), V_LAST(ic), MSP_LAST(ic));
end
fprintf(FileId,'\n');

fprintf(FileId, 'Quantities averaged over time: %d h to %d h:\n', Tstart, Tend);
fprintf(FileId,'\n');
fprintf(FileId, '%15s %10s %10s %10s %10s %10s %10s %10s\n', 'Case', 'IKE', 'RMW', 'R64KT', 'R50KT', 'R34KT', 'V10', 'MSP');
fprintf(FileId,'\n');
for ic = 1:length(CASES)
  fprintf(FileId,'%15s %10.2e %10.2e %10.2e %10.2e %10.2e %10.2e %10.2e\n', ...
    CASES{ic}, KE_AVG(ic), RMW_AVG(ic), R64KT_AVG(ic), R50KT_AVG(ic), R34KT_AVG(ic), V_AVG(ic), MSP_AVG(ic));
end
fprintf(FileId,'\n');

% Close the file if not stdout nor stderr
if (FileId > 2)
  fclose(FileId);
end

end
