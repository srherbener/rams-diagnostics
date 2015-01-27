function [ ] = ReadTogaCoareSound(SoundDir, SoundDate)
% ReadTogaCoareSound read in the TOGA COARE sounding data from version 2 ascii files

  % Data is organized in the ascii files:
  %
  %   <Header>
  %   <One level data>
  %   <Header>
  %   <One level data>
  %   ...
  %
  %  <Header> examples: "T(C) at surface", "T(C) at 1000.mb", "T(C) at 975.mb"
  %
  %  <One level data>:
  %     41 rows     (lon values)
  %     21 columns  (lat values)
  %
  % There are 41 levels (pairs of <Header> and <One level data>).
  %
  % Longitude values (1 degree spacing):
  %    lon(i) = 140 + (i-1)*1, where i = 1:41
  %    lon(i) = 139 + i
  %
  % Latitude values (1 degree spacing):
  %    lat(j) = -10 + (j-1)*1,  where j = 1:21
  %    lat(j) = -11 + j
  %
  % Pressure levels (25 mb spacing):
  %    plev(1) = surface (1013 mb)
  %    plev(k) = 1000 - (k-2)*25, where k = 2:41

  Nlon = 41;
  Nlat = 21;
  Nlev = 41;

  LON = 140:1:180;
  LAT = -10:1:10;
  P = [ 1013 1000:-25:25 ];

  fprintf('Reading TOGA COARE sounding:\n');
  fprintf('  Sounding directory: %s\n', SoundDir);
  fprintf('  Sounding date: %s\n', SoundDate);
  fprintf('\n');

  % read in the pressure file, it contains just the sea level pressure at all lat/lon points (2D)
  InFile = sprintf('%s/p.%s', SoundDir, SoundDate);
  SLP = ReadSoundFile(InFile, Nlon, Nlat, 1);

  % read in T, Q, U and V (3D)
  InFile = sprintf('%s/t.%s', SoundDir, SoundDate);
  T = ReadSoundFile(InFile, Nlon, Nlat, Nlev);

  InFile = sprintf('%s/q.%s', SoundDir, SoundDate);
  Q = ReadSoundFile(InFile, Nlon, Nlat, Nlev);

  InFile = sprintf('%s/u.%s', SoundDir, SoundDate);
  U = ReadSoundFile(InFile, Nlon, Nlat, Nlev);

  InFile = sprintf('%s/v.%s', SoundDir, SoundDate);
  V = ReadSoundFile(InFile, Nlon, Nlat, Nlev);

  fprintf('\n'); 

  % Take horizontal averages of SLP, T, Q, U, V.
  % Data are organized as one of: (x,y) or (x,y,z)
  SLP_AVG = squeeze(mean(mean(SLP,1),2));
  T_AVG   = squeeze(mean(mean(T,1),2));
  Q_AVG   = squeeze(mean(mean(Q,1),2));
  U_AVG   = squeeze(mean(mean(U,1),2));
  V_AVG   = squeeze(mean(mean(V,1),2));

  % Replace P(1) with the SLP average
  P(1) = SLP_AVG;

  % Calculate potential temperature
  %  THETA = T * (P0 / P) ^ (Rd/Cp)
  %
  %    P0 = 1000 mb
  %    Rd = 1004 J/kg/K
  %    Cp = 287  J/kg/K
  %
  %    --> Rd/Cp = 0.286
  %
  %  T needs to be in Kelvin (T + 273.15)
  %
  THETA_AVG = (T_AVG + 273.15) .* ((1000 ./ P') .^ (0.286));

  % Adust the theta profile so that the first entry is 300K.
  %
  %   Create deleta: TH(1) - 300
  %   For first five entries: New TH = TH - delta
  %   Next five entries: New TH = TH - (delta * (10-i)/5)
  %   Remaining entries: New TH = TH
  DELTA_TH = THETA_AVG(1) - 300;
  TH_ADJUST = vertcat( ones([ 5 1]), ([4:-1:0]' ./ 5), zeros([ 31 1 ]) ) .* DELTA_TH;
  THETA_300 = THETA_AVG - TH_ADJUST;

  % Save both the original and average values
  OutFile = 'DIAGS/TogaCoareSoundings.h5';
  fprintf('  Writing: %s\n', OutFile);

  hdf5write(OutFile, 'LON', LON);
  hdf5write(OutFile, 'LAT', LAT, 'WriteMode', 'append');
  hdf5write(OutFile, 'P',   P,   'WriteMode', 'append');
   
  hdf5write(OutFile, 'SLP_ALL', SLP, 'WriteMode', 'append');
  hdf5write(OutFile, 'T_ALL',   T,   'WriteMode', 'append');
  hdf5write(OutFile, 'Q_ALL',   Q,   'WriteMode', 'append');
  hdf5write(OutFile, 'U_ALL',   U,   'WriteMode', 'append');
  hdf5write(OutFile, 'V_ALL',   V,   'WriteMode', 'append');

  hdf5write(OutFile, 'SLP_AVG',   SLP_AVG,   'WriteMode', 'append');
  hdf5write(OutFile, 'T_AVG',     T_AVG,     'WriteMode', 'append');
  hdf5write(OutFile, 'THETA_AVG', THETA_AVG, 'WriteMode', 'append');
  hdf5write(OutFile, 'THETA_300', THETA_300, 'WriteMode', 'append');
  hdf5write(OutFile, 'Q_AVG',     Q_AVG,     'WriteMode', 'append');
  hdf5write(OutFile, 'U_AVG',     U_AVG,     'WriteMode', 'append');
  hdf5write(OutFile, 'V_AVG',     V_AVG,     'WriteMode', 'append');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ReadSoundFile()
%
% This routine will read the contents of a toga coare soudinging file
% into the output, SOUT. SOUT will be a 3D array, (Nlon, Nlat, Nlev).
%
function [ SOUT ] = ReadSoundFile(InFile, Nlon, Nlat, Nlev);

  % See comments above for file organization

  fprintf('  Reading: %s\n', InFile);

  Fid = fopen(InFile);

  % Create a format string that contains Nlat number of %f fields (ie,
  % the number of columns in the tabular data within the file).
  TsFormat = '';
  for i = 1:Nlat
    TsFormat = sprintf('%s %%f', TsFormat); 
  end

  % Read in the data
  %
  % The textscan command below will skip over <Header> due to the
  % 'HeaderLines', 1 spec, then it will read lines matching TsFormat
  % up to the next <Header> line. The result in RAW will a cell array
  % where each cell is a column vector containg the values in the 
  % corresponding column of the tabular format in the input file.
  % This fits nicely with cell2mat which will reconstruct the 2D
  % array for the current level.
  SOUT = zeros([ Nlon Nlat Nlev ]);
  Fid = fopen(InFile);
  for k = 1:Nlev
    RAW = textscan(Fid, TsFormat, 'HeaderLines', 1);
    SOUT(:,:,k) = cell2mat(RAW);
  end

  fclose(Fid);
end

