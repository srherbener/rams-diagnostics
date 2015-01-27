function [ ] = TogaCoareToRamsSound(SoundFile)
% TogaCoareToRamsSound convert TOGA COARE sounding to RAMS format

  fprintf('Converting TOGA COARE sounding:\n');
  fprintf('  Sounding file: %s\n', SoundFile);
  fprintf('\n');

  % Read in the averaged profiles for: P, T, Q, U and V
  % RAMS uses mixing ratio (R), whereas TOGA COARE uses
  % specific humidity. R = Q / (1 - Q) where R and Q are
  % in kg/kg. Q is stored in g/kg.
  fprintf('  Pressure\n');
  PS = h5read(SoundFile, '/P');

  fprintf('  Theta\n');
  TS = h5read(SoundFile, '/THETA_300'); % read theta that was adjusted to 300K surface temp

  fprintf('  Specific humidity\n');
  Q = h5read(SoundFile, '/Q_AVG') ./ 1000;
  RTS = (Q ./ (1 - Q)) .* 1000;

  fprintf('  Zonal wind\n');
  US = h5read(SoundFile, '/U_AVG');

  fprintf('  Meridional wind\n');
  VS = h5read(SoundFile, '/V_AVG');

  fprintf('\n');

  % dump sounding data in RAMS format into the output file
  OutFile = 'DIAGS/TogaCoareRamsSounding.txt';
  fprintf('  Writing: %s\n', OutFile);
  fprintf('\n');

  Fid = fopen(OutFile, 'w');

  DumpRamsSounding(Fid, PS, 'PS');
  DumpRamsSounding(Fid, TS, 'TS');
  DumpRamsSounding(Fid, RTS, 'RTS');
  DumpRamsSounding(Fid, US, 'US');
  DumpRamsSounding(Fid, VS, 'VS');

  fclose(Fid);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DumpRamsSounding()
%
% Write (to stdout) out the sounding in Rams format
%
function [] = DumpRamsSounding(Fid, Prof, Pname)

  Nprof = length(Prof);

  for i = 1:Nprof
    if (i == 1)
      fprintf(Fid, '   %-5s =', Pname);
    end

    fprintf(Fid, ' %10.1f,', Prof(i));

    if (mod(i, 5) == 0)
      fprintf(Fid, '\n');
      fprintf(Fid, '%10s', ' ');
    end
  end

  fprintf(Fid, '\n\n');
end
