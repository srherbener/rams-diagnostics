function [ ] = GenWvarTavg(ConfigFile)
% GenWvarTavg generate text file containing W variance profile

    % Read the config file to get the structure of how the data is laid out in
    % the file system.
    [ Config ] = ReadConfig(ConfigFile);

    Ddir = Config.DiagDir;

    T1 = 145; % 12 h
    T2 = 433; % 36 h

    InVarList  = {
      { 'moments' '/w-w_all'          'avg_wvar'         }
      { 'moments' '/w-w_up0p10_all'   'avg_wvar_up0p10'  }
      { 'moments' '/w-w_dn0p10_all'   'avg_wvar_dn0p10'  }
      { 'moments' '/w-w-w_all'        'avg_wskew'        }
      { 'moments' '/w-w-w_up0p10_all' 'avg_wskew_up0p10' }
      { 'moments' '/w-w-w_dn0p10_all' 'avg_wskew_dn0p10' }
      };

    for icase = 1:length(Config.Cases)
      Case = Config.Cases(icase).Cname;
      fprintf('***************************************************************\n');
      fprintf('Generating time averaged W variance/skew profiles:\n');
      fprintf('  Case: %s\n', Case);
      fprintf('\n');

      for ivar = 1:length(InVarList)
        InFprefix  = InVarList{ivar}{1};
        InVar      = InVarList{ivar}{2};
        OutFprefix = InVarList{ivar}{3};

        InFile = sprintf('%s/%s_%s.h5', Ddir, InFprefix, Case);
        OutFile = sprintf('%s/%s_%s.txt', Ddir, OutFprefix, Case);

        fprintf('  Input file: %s (%s)\n', InFile, InVar);

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
    end
