function [ ] = GenTkeTavg(ConfigFile)
% GenTkeTavg generate time average of Tke profiles

    % Read the config file to get the structure of how the data is laid out in
    % the file system.
    [ Config ] = ReadConfig(ConfigFile);

    Ddir = Config.DiagDir;
    Hdir = 'HDF5';

    T1 = 145; % 12 h
    T2 = 433; % 36 h

    InVar = 'tke_grid';

    for icase = 1:length(Config.Cases)
      Case = Config.Cases(icase).Cname;
      InFile = sprintf('%s/tke-%s-AS-1999-02-10-040000-g1.h5', Hdir, Case);
      OutFile = sprintf('%s/avg_tke_%s.txt', Ddir, Case);

      fprintf('***************************************************************\n');
      fprintf('Generating time averaged TKE profile:\n');
      fprintf('  Case: %s\n', Case);
      fprintf('  Input file: %s (%s)\n', InFile, InVar);
      fprintf('\n');

      % TKE profile is in the i,j = 1,1 column.
      % TKE will be organized as: (t,z,y,x)
      % Z will be organized as: (z)
      TKE_DS = ncgeodataset(InFile);

      TKE_VAR = TKE_DS.geovariable(InVar);
      Z_VAR   = TKE_DS.geovariable('z_coords');

      TKE = squeeze(TKE_VAR.data(T1:T2,:,1,1));
      Z   = squeeze(Z_VAR.data(:));

      % TKE will be organized as: (t,z)
      % averaged over time (dim 1)
      TKE_AVG = squeeze(mean(TKE,1));
      Nz = length(TKE_AVG);

      % write out z and tke values in text form
      fprintf('  Writing: %s\n', OutFile);

      fid = fopen(OutFile, 'wt');
      for i = 1:Nz
        fprintf(fid, '%e %e\n', Z(i), TKE_AVG(i));
      end
      fclose(fid);

      fprintf('\n');
    end
