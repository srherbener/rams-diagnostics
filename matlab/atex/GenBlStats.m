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
  % Third string denotes the type of the remaining input files:
  %     'turb_cov'    --> tsavg, 'turb_cov' function
  %     'gen_moments' --> gen_moments run with 2 input variables
  VarList = {
    { 'cloud_M1_c0p01' 'cloud' 'turb_cov' 'turb_cov_w_theta' 'turb_cov_theta' 'turb_cov_w_vapor' 'turb_cov_vapor' 'turb_cov_w_theta_v' 'turb_cov_theta_v' 'turb_all' }
    { 'cloud_M1_c0p01' 'cloud' 'gen_moments' 'w_theta_flux' 'w-theta' 'w_vapor_flux' 'w-vapor' 'w_theta_v_flux' 'w-theta_v' 'gm_all'}
    { 'cloud_M1_c0p01' 'cloud' 'gen_moments' 'w_theta_flux_stall' 'w-theta' 'w_vapor_flux_stall' 'w-vapor' 'w_theta_v_flux_stall' 'w-theta_v' 'gm_stall' }
    { 'cloud_M1_c0p01' 'cloud' 'gen_moments' 'w_theta_flux_all_cld' 'w-theta' 'w_vapor_flux_all_cld' 'w-vapor' 'w_theta_v_flux_all_cld' 'w-theta_v' 'gm_all_cld' }

% THESE CAUSE problems with theta and consequently Zi, and all We calculations
%  Not sure why: but theta comes out serverly distorted - probably due to the data selection in these being done at discrete points
%  instead of whole columns.
%    { 'cloud_M1_c0p01' 'cloud' 'gen_moments' 'w_theta_flux_ud0p10' 'w-theta' 'w_vapor_flux_ud0p10' 'w-vapor' 'w_theta_v_flux_ud0p10' 'w-theta_v' }
%    { 'cloud_M1_c0p01' 'cloud' 'gen_moments' 'w_theta_flux_up0p10' 'w-theta' 'w_vapor_flux_up0p10' 'w-vapor' 'w_theta_v_flux_up0p10' 'w-theta_v' }
%    { 'cloud_M1_c0p01' 'cloud' 'gen_moments' 'w_theta_flux_dn0p10' 'w-theta' 'w_vapor_flux_dn0p10' 'w-vapor' 'w_theta_v_flux_dn0p10' 'w-theta_v' }

    };

  OutFprefixList = 'bl_stats_0p01';

  NumPtsVar = 'num_points';

  Nfilter = 5;

  NhorizPts = 398 * 398;  % 400 x 400 domain with the outer edges stripped off

  for icase = 1:length(Config.Cases)
    Case = Config.Cases(icase).Cname;
    fprintf('***************************************************************\n');
    fprintf('Generating BL stats:\n');
    fprintf('  Case: %s\n', Case);
  
    OutFile = sprintf('%s/%s_%s.h5', Ddir, OutFprefixList, Case);
    fprintf('  Output file: %s\n', OutFile);
    fprintf('\n');
    hdf5write(OutFile, 'header', 'ATEX');

    for ivar = 1:length(VarList)
      % Always get cloud file from Ddir
      CldFile = sprintf('%s/%s_%s.h5', Ddir, VarList{ivar}{1}, Case);
      CldVar   = VarList{ivar}{2};

      % If the w-theta file has a 'turb_cov_' prefix, assume all the covariance files
      % (remainder on the list) are the tsavg 'turb_cov' type. Otherwise assume they
      % are the gen_moments type generated with 2 variables
      IsTurbType = strcmp(VarList{ivar}{3}, 'turb_cov');
      if(IsTurbType)
        InDir = Tdir;
      else
        InDir = Ddir;
      end

      W_ThFile  = sprintf('%s/%s_%s.h5', InDir, VarList{ivar}{4}, Case);
      W_VapFile = sprintf('%s/%s_%s.h5', InDir, VarList{ivar}{6}, Case);
      W_ThvFile = sprintf('%s/%s_%s.h5', InDir, VarList{ivar}{8}, Case);

      W_ThVar  = VarList{ivar}{5};
      W_VapVar = VarList{ivar}{7};
      W_ThvVar = VarList{ivar}{9};

      VarTag = VarList{ivar}{10};
      
      fprintf('  Input cloud mixing ratio file: %s\n', CldFile);
      fprintf('    Variable name: %s\n', CldVar);
      fprintf('  Input w-theta file: %s\n', W_ThFile);
      fprintf('    Variable name: %s\n', W_ThVar);
      fprintf('  Input w-vapor file: %s\n', W_VapFile);
      fprintf('    Variable name: %s\n', W_VapVar);
      fprintf('  Input w-theta_v file: %s\n', W_ThvFile);
      fprintf('    Variable name: %s\n', W_ThvVar);
      
      % coordinates
      Z = hdf5read(CldFile, 'z_coords');
      T = hdf5read(CldFile, 't_coords');

      Nz = length(Z);
      Nt = length(T);

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
  
      % need to generate the moments differently depending upon the input type
      if (IsTurbType)
        % output from tsavg 'turb_cov'
        W_TH_TERMS = squeeze(hdf5read(W_ThFile, W_ThVar));
        TH_SUMS   = squeeze(W_TH_TERMS(2, :, :));
        W_TH_SUMS = squeeze(W_TH_TERMS(3, :, :));
        W_TH_NPTS = squeeze(W_TH_TERMS(4, :, :));
        TH = TH_SUMS ./ W_TH_NPTS;
        W_TH = W_TH_SUMS ./ W_TH_NPTS;

        W_VAP_TERMS = squeeze(hdf5read(W_VapFile, W_VapVar));
        VAP_SUMS   = squeeze(W_VAP_TERMS(2, :, :));
        W_VAP_SUMS = squeeze(W_VAP_TERMS(3, :, :));
        W_VAP_NPTS = squeeze(W_VAP_TERMS(4, :, :));
        VAP = VAP_SUMS ./ W_VAP_NPTS;
        W_VAP = W_VAP_SUMS ./ W_VAP_NPTS;

        W_THV_TERMS = squeeze(hdf5read(W_ThvFile, W_ThvVar));
        THV_SUMS   = squeeze(W_THV_TERMS(2, :, :));
        W_THV_SUMS = squeeze(W_THV_TERMS(3, :, :));
        W_THV_NPTS = squeeze(W_THV_TERMS(4, :, :));
        THV = THV_SUMS ./ W_THV_NPTS;
        W_THV = W_THV_SUMS ./ W_THV_NPTS;
 
        % need to add code for these - average over time interval
        TH_TALL   = nan([Nz 1]);
        W_TH_TALL = nan([Nz 1]);

        VAP_TALL   = nan([Nz 1]);
        W_VAP_TALL = nan([Nz 1]);

        THV_TALL   = nan([Nz 1]);
        W_THV_TALL = nan([Nz 1]);
      else
        % output from gen_moments using 2 variables
        W_TH_TERMS = squeeze(hdf5read(W_ThFile, W_ThVar));
        W_TH_NPTS  = squeeze(hdf5read(W_ThFile, 'num_points'));
        TH   = zeros([ Nz Nt ]);
        W_TH = zeros([ Nz Nt ]);

        W_VAP_TERMS = squeeze(hdf5read(W_VapFile, W_VapVar));
        W_VAP_NPTS  = squeeze(hdf5read(W_VapFile, 'num_points'));
        VAP   = zeros([ Nz Nt ]);
        W_VAP = zeros([ Nz Nt ]);

        W_THV_TERMS = squeeze(hdf5read(W_ThvFile, W_ThvVar));
        W_THV_NPTS  = squeeze(hdf5read(W_ThvFile, 'num_points'));
        THV   = zeros([ Nz Nt ]);
        W_THV = zeros([ Nz Nt ]);
        
        for i = 1:Nt
          % last argument says whether or not to do averaging across time first (before spatial averaging)
          % since doing one time step per iteration, it's more efficient to set this to '1'
          % MOMENTS will be formed as: (Term, Order, z)
          %    (Term is term number, Order is order of terms)
          %    For two variables: w, theta
          %
          %    Term,  Order-->    1                 2
          %      1              mean(w)        cov(w,theta)
          %      2              mean(theta)    0
          %
          [ MOMENTS OUT_NPTS ] = GenMoments(W_TH_TERMS, W_TH_NPTS, i, i, 1);
          TH(:,i) = squeeze(MOMENTS(2, 1, :));
          W_TH(:,i) = squeeze(MOMENTS(1, 2, :));

          [ MOMENTS OUT_NPTS ] = GenMoments(W_VAP_TERMS, W_VAP_NPTS, i, i, 1);
          VAP(:,i)   = squeeze(MOMENTS(2, 1, :));
          W_VAP(:,i) = squeeze(MOMENTS(1, 2, :));

          [ MOMENTS OUT_NPTS ] = GenMoments(W_THV_TERMS, W_THV_NPTS, i, i, 1);
          THV(:,i)   = squeeze(MOMENTS(2, 1, :));
          W_THV(:,i) = squeeze(MOMENTS(1, 2, :));
        end

        % Get averages over the 24-hour analysis period
        [ MOMENTS OUT_NPTS ] = GenMoments(W_TH_TERMS, W_TH_NPTS, 145, 433, 1);
        TH_TALL = squeeze(MOMENTS(2, 1, :));
        W_TH_TALL = squeeze(MOMENTS(1, 2, :));

        [ MOMENTS OUT_NPTS ] = GenMoments(W_VAP_TERMS, W_VAP_NPTS, 145, 433, 1);
        VAP_TALL   = squeeze(MOMENTS(2, 1, :));
        W_VAP_TALL = squeeze(MOMENTS(1, 2, :));

        [ MOMENTS OUT_NPTS ] = GenMoments(W_THV_TERMS, W_THV_NPTS, 145, 433, 1);
        THV_TALL   = squeeze(MOMENTS(2, 1, :));
        W_THV_TALL = squeeze(MOMENTS(1, 2, :));
      end

      % Profile arrays are organized as: (z,t)
      %
      % For the data that is conditionally sampled at points, it is possible
      % to get zero counts (with an associated zero sum) that will have
      % produced a nan in these arrays. It won't work to just leave the nan
      % in the array, and it won't work to substitute all nans in a given
      % array with a constant value. Instead, use an average of the
      % points above and below the nan.
      TH = InpaintNanColumns(TH);
      W_TH = InpaintNanColumns(W_TH);

      VAP = InpaintNanColumns(VAP);
      W_VAP = InpaintNanColumns(W_VAP);

      THV = InpaintNanColumns(THV);
      W_THV = InpaintNanColumns(W_THV);

      % Find inversion level (Zi) by locating max dtheta/dz in vertical
      %
      [ Zinv ZinvInd ] = FindInversion(TH, Z);

      [ Zinv_TALL ZinvInd_TALL ] = FindInversion(TH, Z);

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
        %   form jump conditions (delta) from ZinvInd-1 and ZinvInd+1 levels
        %   use average of flux from ZinvInd-1 and ZinvInd+1 levels
        %   entrainment velocity is flux / (jump conditions)
        Z1 = ZinvInd(i) - 1;
        if (Z1 < 1)
          Z1 = 1;
        end
        Z2 = ZinvInd(i) + 1;

        % Heat entrainment velocity (theta)
        TH_I = squeeze(TH(:,i));
        W_TH_I = squeeze(W_TH(:,i));
        D_TH = TH_I(Z2) - TH_I(Z1);
        W_TH_FLUX = (W_TH_I(Z2) + W_TH_I(Z1)) * 0.5;
        ThetaWe(i) = W_TH_FLUX / D_TH;

        % Moisture entrainment velocity (vapor)
        VAP_I = squeeze(VAP(:,i));
        W_VAP_I = squeeze(W_VAP(:,i));
        D_VAP = VAP_I(Z2) - VAP_I(Z1);
        W_VAP_FLUX = (W_VAP_I(Z2) + W_VAP_I(Z1)) * 0.5;
        VaporWe(i) = W_VAP_FLUX / D_VAP;

        % Bouyancy entrainment velocity (theta_v)
        THV_I = squeeze(THV(:,i));
        W_THV_I = squeeze(W_THV(:,i));
        D_THV = THV_I(Z2) - THV_I(Z1);
        W_THV_FLUX = (W_THV_I(Z2) + W_THV_I(Z1)) * 0.5;
        ThetaV_We(i) = W_THV_FLUX / D_THV;

        if (i == Nt)
          % Time averaged values
          Z1_TALL = ZinvInd_TALL(i) - 1;
          if (Z1_TALL < 1)
            Z1_TALL = 1;
          end
          Z2_TALL = ZinvInd_TALL(i) + 1;

          % Heat entrainment velocity (theta)
          D_TH_TALL = TH_TALL(Z2_TALL) - TH_TALL(Z1_TALL);
          W_TH_FLUX_TALL= (W_TH_TALL(Z2_TALL) + W_TH_TALL(Z1_TALL)) * 0.5;
          ThetaWe_TALL = W_TH_FLUX_TALL / D_TH_TALL;

          % Moisture entrainment velocity (vapor)
          D_VAP_TALL = VAP_TALL(Z2_TALL) - VAP_TALL(Z1_TALL);
          W_VAP_FLUX_TALL= (W_VAP_TALL(Z2_TALL) + W_VAP_TALL(Z1_TALL)) * 0.5;
          VaporWe_TALL = W_VAP_FLUX_TALL / D_VAP_TALL;

          % Bouyancy entrainment velocity (theta_v)
          D_THV_TALL = THV_TALL(Z2_TALL) - THV_TALL(Z1_TALL);
          W_THV_FLUX_TALL= (W_THV_TALL(Z2_TALL) + W_THV_TALL(Z1_TALL)) * 0.5;
          ThetaV_We_TALL = W_THV_FLUX_TALL / D_THV_TALL;
        end

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
      OutVar = reshape(Zinv, [ 1 1 1 Nt ]);
      OutVname = sprintf('InversionHeight_%s', VarTag);
      hdf5write(OutFile, OutVname, OutVar, 'WriteMode', 'append');
      OutVar = reshape(ZinvSmooth, [ 1 1 1 Nt ]);
      OutVname = sprintf('InversionHeightSmooth_%s', VarTag);
      hdf5write(OutFile, OutVname, OutVar, 'WriteMode', 'append');

      OutVar = reshape(ZtopCf, [ 1 1 1 Nt ]);
      OutVname = sprintf('CloudTopCf_%s', VarTag);
      hdf5write(OutFile, OutVname, OutVar, 'WriteMode', 'append');
      OutVar = reshape(ZtopCfSmooth, [ 1 1 1 Nt ]);
      OutVname = sprintf('CloudTopCfSmooth_%s', VarTag);
      hdf5write(OutFile, OutVname, OutVar, 'WriteMode', 'append');

      OutVar = reshape(ZbotCf, [ 1 1 1 Nt ]);
      OutVname = sprintf('CloudBotCf_%s', VarTag);
      hdf5write(OutFile, OutVname, OutVar, 'WriteMode', 'append');
      OutVar = reshape(ZbotCfSmooth, [ 1 1 1 Nt ]);
      OutVname = sprintf('CloudBotCfSmooth_%s', VarTag);
      hdf5write(OutFile, OutVname, OutVar, 'WriteMode', 'append');

      OutVar = reshape(ZtopCld, [ 1 1 1 Nt ]);
      OutVname = sprintf('CloudTopCld_%s', VarTag);
      hdf5write(OutFile, OutVname, OutVar, 'WriteMode', 'append');
      OutVar = reshape(ZtopCldSmooth, [ 1 1 1 Nt ]);
      OutVname = sprintf('CloudTopCldSmooth_%s', VarTag);
      hdf5write(OutFile, OutVname, OutVar, 'WriteMode', 'append');

      OutVar = reshape(ZbotCld, [ 1 1 1 Nt ]);
      OutVname = sprintf('CloudBotCld_%s', VarTag);
      hdf5write(OutFile, OutVname, OutVar, 'WriteMode', 'append');
      OutVar = reshape(ZbotCldSmooth, [ 1 1 1 Nt ]);
      OutVname = sprintf('CloudBotCldSmooth_%s', VarTag);
      hdf5write(OutFile, OutVname, OutVar, 'WriteMode', 'append');

      OutVar = reshape(ThetaWe, [ 1 1 1 Nt ]);
      OutVname = sprintf('ThetaWe_%s', VarTag);
      hdf5write(OutFile, OutVname, OutVar, 'WriteMode', 'append');
      OutVar = reshape(ThetaWeSmooth, [ 1 1 1 Nt ]);
      OutVname = sprintf('ThetaWeSmooth_%s', VarTag);
      hdf5write(OutFile, OutVname, OutVar, 'WriteMode', 'append');

      OutVar = reshape(VaporWe, [ 1 1 1 Nt ]);
      OutVname = sprintf('VaporWe_%s', VarTag);
      hdf5write(OutFile, OutVname, OutVar, 'WriteMode', 'append');
      OutVar = reshape(VaporWeSmooth, [ 1 1 1 Nt ]);
      OutVname = sprintf('VaporWeSmooth_%s', VarTag);
      hdf5write(OutFile, OutVname, OutVar, 'WriteMode', 'append');

      OutVar = reshape(ThetaV_We, [ 1 1 1 Nt ]);
      OutVname = sprintf('ThetaV_We_%s', VarTag);
      hdf5write(OutFile, OutVname, OutVar, 'WriteMode', 'append');
      OutVar = reshape(ThetaV_WeSmooth, [ 1 1 1 Nt ]);
      OutVname = sprintf('ThetaV_WeSmooth_%s', VarTag);
      hdf5write(OutFile, OutVname, OutVar, 'WriteMode', 'append');

      OutVar = reshape(CfMax, [ 1 1 1 Nt ]);
      OutVname = sprintf('CfMax_%s', VarTag);
      hdf5write(OutFile, OutVname, OutVar, 'WriteMode', 'append');
      OutVar = reshape(CfMaxSmooth, [ 1 1 1 Nt ]);
      OutVname = sprintf('CfMaxSmooth_%s', VarTag);
      hdf5write(OutFile, OutVname, OutVar, 'WriteMode', 'append');

      OutVar = reshape(CfMaxHeight, [ 1 1 1 Nt ]);
      OutVname = sprintf('CfMaxHeight_%s', VarTag);
      hdf5write(OutFile, OutVname, OutVar, 'WriteMode', 'append');
      OutVar = reshape(CfMaxHeightSmooth, [ 1 1 1 Nt ]);
      OutVname = sprintf('CfMaxHeightSmooth_%s', VarTag);
      hdf5write(OutFile, OutVname, OutVar, 'WriteMode', 'append');

      TallVarTag = sprintf('%s_TALL', VarTag);
      OutVname = sprintf('ThetaWe_%s', TallVarTag);
      hdf5write(OutFile, OutVname, ThetaWe_TALL, 'WriteMode', 'append');
      OutVname = sprintf('VaporWe_%s', TallVarTag);
      hdf5write(OutFile, OutVname, VaporWe_TALL, 'WriteMode', 'append');
      OutVname = sprintf('ThetaV_We_%s', TallVarTag);
      hdf5write(OutFile, OutVname, ThetaV_We_TALL, 'WriteMode', 'append');

      OutVname = sprintf('TH_%s', TallVarTag);
      hdf5write(OutFile, OutVname, TH_TALL, 'WriteMode', 'append');
      OutVname = sprintf('VAP_%s', TallVarTag);
      hdf5write(OutFile, OutVname, VAP_TALL, 'WriteMode', 'append');
      OutVname = sprintf('THV_%s', TallVarTag);
      hdf5write(OutFile, OutVname, THV_TALL, 'WriteMode', 'append');

      OutVname = sprintf('TH_%s', VarTag);
      hdf5write(OutFile, OutVname, TH, 'WriteMode', 'append');
      OutVname = sprintf('VAP_%s', VarTag);
      hdf5write(OutFile, OutVname, VAP, 'WriteMode', 'append');
      OutVname = sprintf('THV_%s', VarTag);
      hdf5write(OutFile, OutVname, THV, 'WriteMode', 'append');

      fprintf('\n');
    end
    hdf5write(OutFile, 'x_coords', Xdummy, 'WriteMode', 'append');
    hdf5write(OutFile, 'y_coords', Ydummy, 'WriteMode', 'append');
    hdf5write(OutFile, 'z_coords', Z,      'WriteMode', 'append');
    hdf5write(OutFile, 't_coords', T,      'WriteMode', 'append');
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%
function [ Zinv ZinvInd ] = FindInversion(TH, Z)
% FindInversion given theta and z, find the region of maximum dtheta/dz

  % TH is (z,t)
  [ Nz Nt ] = size(TH);

  % Remove kinks from theta so the subsequent find does not get confused
  TH_SMOOTH = zeros([ Nz Nt ]);
  for i = 1:Nt
    TH_SMOOTH(:,i) = smooth(TH(:,i), 3);
  end

  DTH = TH_SMOOTH(2:end, :) - TH_SMOOTH(1:end-1, :);
  DZ  = Z(2:end) - Z(1:end-1);
  DZ_FULL = repmat(DZ, [ 1 Nt ]);
  DTH_DZ = DTH ./ DZ_FULL;

  % ZinvInd will be index of Z at the bottom of the interval where
  % max DthDz exists.
  ZinvInd = zeros([ 1 Nt ]);
  Zinv    = zeros([ 1 Nt ]);
  for i = 1:Nt
    [ Val ZinvInd(i) ] = max(DTH_DZ(:,i));
    Zinv(i) = Z(ZinvInd(i));
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%
function [ OUT ] = InpaintNanColumns(IN)
% InpaintNanColumns call inpaintn on each column separately in a 2d array

% inpaintn is a function that replaces nans with values interpolated from it's neighbors

  [ Nr Nc ] = size(IN);

  OUT = zeros([ Nr Nc ]);
  for i = 1:Nc
    COL = squeeze(IN(:,i));
    if (sum(~isnan(COL)) == 0)
      % all entries are nan, just set to zeros
      OUT(:,i) = zeros([ Nr 1 ]);
    else
      OUT(:,i) = inpaintn(COL,10);
    end
  end
end
