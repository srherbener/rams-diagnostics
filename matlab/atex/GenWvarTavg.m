function [ ] = GenWvarTavg(ConfigFile)
% GenWvarTavg generate text file containing W variance profile

    % Read the config file to get the structure of how the data is laid out in
    % the file system.
    [ Config ] = ReadConfig(ConfigFile);

    Ddir = Config.DiagDir;

    T1 = 145; % 12 h
    T2 = 433; % 36 h

    InVar = '/w-w_all';

    for icase = 1:length(Config.Cases)
      Case = Config.Cases(icase).Cname;
      InFile = sprintf('%s/moments_%s.h5', Ddir, Case);
      OutFile = sprintf('%s/avg_wvar_%s.txt', Ddir, Case);

      fprintf('***************************************************************\n');
      fprintf('Generating time averaged W variance profile:\n');
      fprintf('  Case: %s\n', Case);
      fprintf('  Input file: %s (%s)\n', InFile, InVar);
      fprintf('\n');

      WVAR = squeeze(h5read(InFile, InVar));
      Z   = squeeze(h5read(InFile, '/z_coords'));
      Nz = length(Z);

      % write out z and tke values in text form
      fprintf('  Writing: %s\n', OutFile);

      fid = fopen(OutFile, 'wt');
      for i = 1:Nz
        fprintf(fid, '%e %e\n', Z(i), WVAR(i));
      end
      fclose(fid);

      fprintf('\n');
    end
