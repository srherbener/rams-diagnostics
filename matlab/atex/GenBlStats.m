function [ ] = GenBlStats(ConfigFile)
% GenBlStats generate measurements of boundary layer - inversion, cloud top, cloud bottom

  % Read the config file to get the structure of how the data is laid out in
  % the file system.
  [ Config ] = ReadConfig(ConfigFile);

  Ddir = Config.DiagDir;
  Tdir = Config.TsavgDir;

  % Generate
  %   No filter on theta
  %   Filter clouds based on cloud mix ratios >= 0.01 g/kg
  VarList = {
    { 'cloud_M1_c0p01' 'cloud' 'turb_cov_w_theta' 'turb_cov_theta' 'turb_cov_w_vapor' 'turb_cov_vapor' 'turb_cov_w_theta_v' 'turb_cov_theta_v' }
    };

  OutFprefixList = {
    'bl_stats_0p01'
    };

  NumPtsVar = 'num_points';

  Nfilter = 5;

  NhorizPts = 398 * 398;  % 400 x 400 domain with the outer edges stripped off

  for icase = 1:length(Config.Cases)
    for ivar = 1:length(VarList)
      Case = Config.Cases(icase).Cname;
  
      CldFile = sprintf('%s/%s_%s.h5', Ddir, VarList{ivar}{1}, Case);
      CldVar  = VarList{ivar}{2};
      W_ThFile  = sprintf('%s/%s_%s.h5', Tdir, VarList{ivar}{3}, Case);
      W_ThVar   = VarList{ivar}{4};
      W_VapFile  = sprintf('%s/%s_%s.h5', Tdir, VarList{ivar}{5}, Case);
      W_VapVar   = VarList{ivar}{6};
      W_ThvFile  = sprintf('%s/%s_%s.h5', Tdir, VarList{ivar}{7}, Case);
      W_ThvVar   = VarList{ivar}{8};
      
      OutFile = sprintf('%s/%s_%s.h5', Ddir, OutFprefixList{ivar}, Case);
  
      fprintf('***************************************************************\n');
      fprintf('Generating BL stats:\n');
      fprintf('  Case: %s\n', Case);
      fprintf('  Input cloud mixing ratio file: %s\n', CldFile);
      fprintf('    Variable name: %s\n', CldVar);
      fprintf('  Input w-theta file: %s\n', W_ThFile);
      fprintf('    Variable name: %s\n', W_ThVar);
      fprintf('  Input w-vapor file: %s\n', W_VapFile);
      fprintf('    Variable name: %s\n', W_VapVar);
      fprintf('  Input w-theta_v file: %s\n', W_ThvFile);
      fprintf('    Variable name: %s\n', W_ThvVar);
      
      fprintf('  Output file: %s\n', OutFile);
      fprintf('\n');
  
      % Cloud, Vapor mixing ratio is in g/kg
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
 
      CF = CLD_NPTS ./ NhorizPts; % cloud fraction per level
  
      W_TH_TERMS = squeeze(hdf5read(W_ThFile, W_ThVar));
      TH_SUMS   = squeeze(W_TH_TERMS(2, :, :));
      W_TH_SUMS = squeeze(W_TH_TERMS(3, :, :));
      W_TH_NPTS = squeeze(W_TH_TERMS(4, :, :));
      TH = TH_SUMS ./ W_TH_NPTS;
      TH(isnan(TH)) = 0;
      W_TH = W_TH_SUMS ./ W_TH_NPTS;
      W_TH(isnan(W_TH)) = 0;

      W_VAP_TERMS = squeeze(hdf5read(W_VapFile, W_VapVar));
      VAP_SUMS   = squeeze(W_VAP_TERMS(2, :, :));
      W_VAP_SUMS = squeeze(W_VAP_TERMS(3, :, :));
      W_VAP_NPTS = squeeze(W_VAP_TERMS(4, :, :));
      VAP = VAP_SUMS ./ W_VAP_NPTS;
      VAP(isnan(VAP)) = 0;
      W_VAP = W_VAP_SUMS ./ W_VAP_NPTS;
      W_VAP(isnan(W_VAP)) = 0;

      W_THV_TERMS = squeeze(hdf5read(W_ThvFile, W_ThvVar));
      THV_SUMS   = squeeze(W_THV_TERMS(2, :, :));
      W_THV_SUMS = squeeze(W_THV_TERMS(3, :, :));
      W_THV_NPTS = squeeze(W_THV_TERMS(4, :, :));
      THV = THV_SUMS ./ W_THV_NPTS;
      THV(isnan(THV)) = 0;
      W_THV = W_THV_SUMS ./ W_THV_NPTS;
      W_THV(isnan(W_THV)) = 0;


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
      ZtopCf = zeros([ 1 Nt ]);
      ZbotCf = zeros([ 1 Nt ]);
      ZtopCld = zeros([ 1 Nt ]);
      ZbotCld = zeros([ 1 Nt ]);
      ThetaWe = zeros([ 1 Nt ]);
      VaporWe = zeros([ 1 Nt ]);
      ThetaV_We = zeros([ 1 Nt ]);
      MaxCf = zeros([ 1 Nt ]);
      MaxCfHeight = zeros([ 1 Nt ]);
      for i = 1:Nt
        % inversion
        % use smooth() to eliminate "blips" in Dtheta/Dz, ie make the largest inversion region stand out
        DTH_DZ_I = smooth(squeeze(DTH_DZ(:,i)), Nfilter);
        DTH_DZ_MAX = max(DTH_DZ_I);
        ZinvInd = find(DTH_DZ_I == DTH_DZ_MAX, 1 , 'first'); % Index is set to the bottom of the interval
                                                             % where max Dtheta/Dz was found

        Zinv(i) = (Z(ZinvInd) + Z(ZinvInd+1)) / 2; % Set Zi to average of bottom and top of interval

        % For cloud top and bottom:
        % Find top and bottom using cloud fraction, and cloud mixing ratio
        % Find the max value of both of these and use 0.5 * Max as threshold to find top and bottom
        % Smooth the profiles so that the larger peaks will be emphasized
        %   
        CF_I = smooth(squeeze(CF(:,i)), Nfilter);
        CLD_I = smooth(squeeze(CLD(:,i)), Nfilter);

        CfThreshold = max(CF_I) * 0.1;
        CldThreshold = 0.1; % max(CLD_I) * 0.5;

        ZbotCfInd = find(CF_I >= CfThreshold, 1, 'first');
        ZtopCfInd = find(CF_I >= CfThreshold, 1, 'last');
        ZbotCldInd = find(CLD_I >= CldThreshold, 1, 'first');
        ZtopCldInd = find(CLD_I >= CldThreshold, 1, 'last');

        if (isempty(ZbotCfInd))
          ZbotCfInd = 1;
        end
        if (isempty(ZtopCfInd))
          ZtopCfInd = 1;
        end
        if (isempty(ZbotCldInd))
          ZbotCldInd = 1;
        end
        if (isempty(ZtopCldInd))
          ZtopCldInd = 1;
        end

        ZtopCf(i) = Z(ZtopCfInd);
        ZbotCf(i) = Z(ZbotCfInd);
        ZtopCld(i) = Z(ZtopCldInd);
        ZbotCld(i) = Z(ZbotCldInd);

        % entrainment velocity
        %   form jump conditions (delta) from ZinvInd and ZinvInd+1 levels
        %   use average of flux from ZinvInd and ZinvInd+1 levels
        %   entrainment velocity is flux / (jump conditions)

        % Heat entrainment velocity (theta)
        TH_I = squeeze(TH(:,i));
        W_TH_I = squeeze(W_TH(:,i));
        D_TH = TH_I(ZinvInd+1) - TH_I(ZinvInd);
        W_TH_FLUX = (W_TH_I(ZinvInd+1) + W_TH_I(ZinvInd)) * 0.5;
        ThetaWe(i) = W_TH_FLUX / D_TH;

        % Moisture entrainment velocity (vapor)
        VAP_I = squeeze(VAP(:,i));
        W_VAP_I = squeeze(W_VAP(:,i));
        D_VAP = VAP_I(ZinvInd+1) - VAP_I(ZinvInd);
        W_VAP_FLUX = (W_VAP_I(ZinvInd+1) + W_VAP_I(ZinvInd)) * 0.5;
        VaporWe(i) = W_VAP_FLUX / D_VAP;

        % Bouyancy entrainment velocity (theta_v)
        THV_I = squeeze(THV(:,i));
        W_THV_I = squeeze(W_THV(:,i));
        D_THV = THV_I(ZinvInd+1) - THV_I(ZinvInd);
        W_THV_FLUX = (W_THV_I(ZinvInd+1) + W_THV_I(ZinvInd)) * 0.5;
        ThetaV_We(i) = W_THV_FLUX / D_THV;

        % Find the max value of cloud fraction, also record the height at which this occurs
        CF_I = squeeze(CF(:,i));
        CfMax(i) = max(CF_I);
        ZcfMaxInd = find(CF_I == CfMax(i), 1, 'first');
        CfMaxHeight(i) = Z(ZcfMaxInd);
      end

      % Create smoothed versions of measurements
      ZinvSmooth = smooth(Zinv, Nfilter);
      ZtopCfSmooth = smooth(ZtopCf, Nfilter);
      ZbotCfSmooth = smooth(ZbotCf, Nfilter);
      ZtopCldSmooth = smooth(ZtopCld, Nfilter);
      ZbotCldSmooth = smooth(ZbotCld, Nfilter);
      ThetaWeSmooth = smooth(ThetaWe, Nfilter);
      VaporWeSmooth = smooth(VaporWe, Nfilter);
      ThetaV_WeSmooth = smooth(ThetaV_We, Nfilter);
      CfMaxSmooth = smooth(CfMax, Nfilter);
      CfMaxHeightSmooth = smooth(CfMaxHeight, Nfilter);

      % Save the diameters (in Y) and characteristic diameters (convert m to um)
      Xdummy = 1;
      Ydummy = 1;
      Zdummy = 1;
      OutVar = reshape(Zinv, [ 1 1 1 Nt ]);
      hdf5write(OutFile, 'InversionHeight', OutVar);
      OutVar = reshape(ZinvSmooth, [ 1 1 1 Nt ]);
      hdf5write(OutFile, 'InversionHeightSmooth', OutVar, 'WriteMode', 'append');

      OutVar = reshape(ZtopCf, [ 1 1 1 Nt ]);
      hdf5write(OutFile, 'CloudTopCf', OutVar, 'WriteMode', 'append');
      OutVar = reshape(ZtopCfSmooth, [ 1 1 1 Nt ]);
      hdf5write(OutFile, 'CloudTopCfSmooth', OutVar, 'WriteMode', 'append');

      OutVar = reshape(ZbotCf, [ 1 1 1 Nt ]);
      hdf5write(OutFile, 'CloudBotCf', OutVar, 'WriteMode', 'append');
      OutVar = reshape(ZbotCfSmooth, [ 1 1 1 Nt ]);
      hdf5write(OutFile, 'CloudBotCfSmooth', OutVar, 'WriteMode', 'append');

      OutVar = reshape(ZtopCld, [ 1 1 1 Nt ]);
      hdf5write(OutFile, 'CloudTopCld', OutVar, 'WriteMode', 'append');
      OutVar = reshape(ZtopCldSmooth, [ 1 1 1 Nt ]);
      hdf5write(OutFile, 'CloudTopCldSmooth', OutVar, 'WriteMode', 'append');

      OutVar = reshape(ZbotCld, [ 1 1 1 Nt ]);
      hdf5write(OutFile, 'CloudBotCld', OutVar, 'WriteMode', 'append');
      OutVar = reshape(ZbotCldSmooth, [ 1 1 1 Nt ]);
      hdf5write(OutFile, 'CloudBotCldSmooth', OutVar, 'WriteMode', 'append');

      OutVar = reshape(ThetaWe, [ 1 1 1 Nt ]);
      hdf5write(OutFile, 'ThetaWe', OutVar, 'WriteMode', 'append');
      OutVar = reshape(ThetaWeSmooth, [ 1 1 1 Nt ]);
      hdf5write(OutFile, 'ThetaWeSmooth', OutVar, 'WriteMode', 'append');

      OutVar = reshape(VaporWe, [ 1 1 1 Nt ]);
      hdf5write(OutFile, 'VaporWe', OutVar, 'WriteMode', 'append');
      OutVar = reshape(VaporWeSmooth, [ 1 1 1 Nt ]);
      hdf5write(OutFile, 'VaporWeSmooth', OutVar, 'WriteMode', 'append');

      OutVar = reshape(ThetaV_We, [ 1 1 1 Nt ]);
      hdf5write(OutFile, 'ThetaV_We', OutVar, 'WriteMode', 'append');
      OutVar = reshape(ThetaV_WeSmooth, [ 1 1 1 Nt ]);
      hdf5write(OutFile, 'ThetaV_WeSmooth', OutVar, 'WriteMode', 'append');

      OutVar = reshape(CfMax, [ 1 1 1 Nt ]);
      hdf5write(OutFile, 'CfMax', OutVar, 'WriteMode', 'append');
      OutVar = reshape(CfMaxSmooth, [ 1 1 1 Nt ]);
      hdf5write(OutFile, 'CfMaxSmooth', OutVar, 'WriteMode', 'append');

      OutVar = reshape(CfMaxHeight, [ 1 1 1 Nt ]);
      hdf5write(OutFile, 'CfMaxHeight', OutVar, 'WriteMode', 'append');
      OutVar = reshape(CfMaxHeightSmooth, [ 1 1 1 Nt ]);
      hdf5write(OutFile, 'CfMaxHeightSmooth', OutVar, 'WriteMode', 'append');

      hdf5write(OutFile, 'x_coords', Xdummy, 'WriteMode', 'append');
      hdf5write(OutFile, 'y_coords', Ydummy, 'WriteMode', 'append');
      hdf5write(OutFile, 'z_coords', Zdummy, 'WriteMode', 'append');
      hdf5write(OutFile, 't_coords', T,      'WriteMode', 'append');
      fprintf('\n');
    end
  end
end
