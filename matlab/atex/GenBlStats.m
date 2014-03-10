function [ ] = GenBlStats(ConfigFile)
% GenBlStats generate measurements of boundary layer - inversion, cloud top, cloud bottom

  % Read the config file to get the structure of how the data is laid out in
  % the file system.
  [ Config ] = ReadConfig(ConfigFile);

  Ddir = Config.DiagDir;

  % Generate two cases:
  %   No filter
  %   Filter based on cloud mix ratios >= 0.01 g/kg
  VarList = {
    { 'cloud_M1'       'cloud' 'w_theta_flux' 'w-theta' 2 1 }
    { 'cloud_M1_c0p01' 'cloud' 'w_theta_flux' 'w-theta' 2 1 }
    };

  OutFprefixList = {
    'bl_stats'
    'bl_stats_0p01'
    };

  NumPtsVar = 'num_points';
  Nfilter = 5;

  for icase = 1:length(Config.Cases)
    for ivar = 1:length(VarList)
      Case = Config.Cases(icase).Cname;
  
      CldFile = sprintf('%s/%s_%s.h5', Ddir, VarList{ivar}{1}, Case);
      CldVar  = VarList{ivar}{2};
      ThFile  = sprintf('%s/%s_%s.h5', Ddir, VarList{ivar}{3}, Case);
      ThVar   = VarList{ivar}{4};
      ThNsum  = VarList{ivar}{5};
      ThOrder = VarList{ivar}{6};
      
      OutFile = sprintf('%s/%s_%s.h5', Ddir, OutFprefixList{ivar}, Case);
  
      fprintf('***************************************************************\n');
      fprintf('Generating BL stats:\n');
      fprintf('  Case: %s\n', Case);
      fprintf('  Input cloud mixing ratio file: %s\n', CldFile);
      fprintf('    Variable name: %s\n', CldVar);
      fprintf('  Input theta file: %s\n', ThFile);
      fprintf('    Variable name: %s\n', ThVar);
      fprintf('    Sum number: %d\n', ThNsum);
      fprintf('    Order: %d\n', ThOrder);
      
      fprintf('  Output file: %s\n', OutFile);
      fprintf('\n');
  
      % Cloud mixing ratio is in g/kg
      % Theta is in K
      %
      % In places where the number of points is zero, the sum will also be zero.
      % Dividing sum by number of points in these cases will do 0/0 which will
      % produce a nan. Change the nans back to zeros since want zero result for
      % zero sums.
      %
      % CLD and TH will be organized as (z,t)
      %
      CLD_SUMS = squeeze(hdf5read(CldFile, CldVar));
      CLD_NPTS = squeeze(hdf5read(CldFile, NumPtsVar));
      CLD = CLD_SUMS ./ CLD_NPTS;
      CLD(isnan(CLD)) = 0;
  
      TH_TERMS = squeeze(hdf5read(ThFile, ThVar));
      TH_NPTS  = squeeze(hdf5read(ThFile, NumPtsVar));
      TH_SUMS  = squeeze(TH_TERMS(ThNsum, ThOrder, :, :));
      TH = TH_SUMS ./ TH_NPTS;
      TH(isnan(TH)) = 0;

      % coordinates
      Z = hdf5read(CldFile, 'z_coords');
      T = hdf5read(CldFile, 't_coords');

      Nz = length(Z);
      Nt = length(T);

      % Find inversion level (Zi) by locating max dtheta/dz in vertical
      %
      DTH = TH(2:end, :) - TH(1:end-1, :);
      DZ  = Z(2:end) - Z(1:end-1);
      DZ_FULL = repmat(DZ, [ 1 Nt ]);
      DTH_DZ = DTH ./ DZ_FULL;

      Zinv = zeros([ 1 Nt ]);
      Ztop = zeros([ 1 Nt ]);
      Zbot = zeros([ 1 Nt ]);
      for i = 1:Nt
        % inversion
        DTH_DZ_I = squeeze(DTH_DZ(:,i));
        DTH_DZ_MAX = max(DTH_DZ_I);
        Z1 = find(DTH_DZ_I == DTH_DZ_MAX, 1 , 'first');

        Zinv(i) = Z(Z1);

        % For cloud top and bottom:
        %   1) Find the max cloud mix ratio in the column
        %   2) Set CldThreshold = to 0.5 CldMax
        %   3) Bottom is first occurrence of mix ratio >= CldThreshold
        %   4) Top is last occurrence of mix ratio >= CldThreshold
        CLD_I = squeeze(CLD(:,i));
        CldThreshold = max(CLD_I) .* 0.5;

        % cloud bottom
        Z1 = find(CLD_I >= CldThreshold, 1, 'first');
        if (isempty(Z1))
          % no clouds in this column
          Zbot(i) = 0;
        else
          Zbot(i) = Z(Z1);
        end

        % cloud top
        Z1 = find(CLD_I >= CldThreshold, 1, 'last');
        if (isempty(Z1))
          % no clouds in this column
          Zbot(i) = 0;
        else
          Ztop(i) = Z(Z1);
        end
      end

      % Create smoothed versions of measurements
      ZinvSmooth = smooth(Zinv, 5);
      ZtopSmooth = smooth(Ztop, 5);
      ZbotSmooth = smooth(Zbot, 5);

      % Save the diameters (in Y) and characteristic diameters (convert m to um)
      Xdummy = 1;
      Ydummy = 1;
      Zdummy = 1;
      OutVar = reshape(Zinv, [ 1 1 1 Nt ]);
      hdf5write(OutFile, 'InversionHeight', OutVar);
      OutVar = reshape(ZinvSmooth, [ 1 1 1 Nt ]);
      hdf5write(OutFile, 'InversionHeightSmooth', OutVar, 'WriteMode', 'append');

      OutVar = reshape(Ztop, [ 1 1 1 Nt ]);
      hdf5write(OutFile, 'CloudTop', OutVar, 'WriteMode', 'append');
      OutVar = reshape(ZtopSmooth, [ 1 1 1 Nt ]);
      hdf5write(OutFile, 'CloudTopSmooth', OutVar, 'WriteMode', 'append');

      OutVar = reshape(Zbot, [ 1 1 1 Nt ]);
      hdf5write(OutFile, 'CloudBot', OutVar, 'WriteMode', 'append');
      OutVar = reshape(ZbotSmooth, [ 1 1 1 Nt ]);
      hdf5write(OutFile, 'CloudBotSmooth', OutVar, 'WriteMode', 'append');

      hdf5write(OutFile, 'x_coords', Xdummy, 'WriteMode', 'append');
      hdf5write(OutFile, 'y_coords', Ydummy, 'WriteMode', 'append');
      hdf5write(OutFile, 'z_coords', Zdummy, 'WriteMode', 'append');
      hdf5write(OutFile, 't_coords', T,      'WriteMode', 'append');
      fprintf('\n');
    end
  end
end
