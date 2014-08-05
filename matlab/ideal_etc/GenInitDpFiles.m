function [ ] = GenInitDpFiles(Case, DpFile)
% GenInitDpFiles generate dp files containing initial conditions for Thorncroft et al., 1993 ETC lifecycles

  % Norhtern Hemisphere winter conditions, pick an arbitrary date in winter
  Year = 2014;
  Month = 01;
  Day = 01;
  Time = 1200;

  fprintf('Generating DP files with ETC initial conditions:\n');
  fprintf('  Case: %d\n', Case);
  fprintf('  Writing: %s\n', DpFile);

  % Generate the fields, then dump into the file
  % Use the 'W' (upper case W) on fopen to avoid "out of memory" issues, also
  % fprintf's run much faster when using 'W'
  Fid = fopen(DpFile, 'Wt');
  [ U V T Zg RH Lon Lat Press ] = GenInitFields(Case);
  WriteDpFile(U, V, T, Zg, RH, Lon, Lat, Press, Year, Month, Day, Time, Fid);
  fclose(Fid);
end
