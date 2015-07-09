function [ ] = ExtractSampleProfiles()
% ExtractSampleProfiles extract vertical profile samples

  Ddir = 'DIAGS';

  % make sure output directory exists
  if (exist(Ddir, 'dir') ~= 7)
    mkdir(Ddir);
  end

  % input_file input_dataset case_name outfile_prefix output_dataset
  VarSets = {
      {
        { 'HDF5/RCE_EXP_S70LY/HDF5/tempc-RCE_EXP_S70LY-AC-2012-01-01-000000-g1.h5' '/tempc' }
        { 'HDF5/RCE_EXP_S70LY/HDF5/theta-RCE_EXP_S70LY-AC-2012-01-01-000000-g1.h5' '/theta' }
        { 'HDF5/RCE_EXP_S70LY/HDF5/cloud-RCE_EXP_S70LY-AC-2012-01-01-000000-g1.h5' '/cloud' }
        { 'HDF5/RCE_EXP_S70LY/HDF5/vint_vapor-RCE_EXP_S70LY-AC-2012-01-01-000000-g1.h5' '/vertint_vapor' }
        'RCE_EXP_S70LY'
      }

      {
        { 'HDF5/RCE_EXP_S70LN/HDF5/tempc-RCE_EXP_S70LN-AC-2012-01-01-000000-g1.h5' '/tempc' }
        { 'HDF5/RCE_EXP_S70LN/HDF5/theta-RCE_EXP_S70LN-AC-2012-01-01-000000-g1.h5' '/theta' }
        { 'HDF5/RCE_EXP_S70LN/HDF5/cloud-RCE_EXP_S70LN-AC-2012-01-01-000000-g1.h5' '/cloud' }
        { 'HDF5/RCE_EXP_S70LN/HDF5/vint_vapor-RCE_EXP_S70LN-AC-2012-01-01-000000-g1.h5' '/vertint_vapor' }
        'RCE_EXP_S70LN'
      }

      {
        { 'HDF5/RCE_EXP_S50LN/HDF5/tempc-RCE_EXP_S50LN-AC-2012-01-01-000000-g1.h5' '/tempc' }
        { 'HDF5/RCE_EXP_S50LN/HDF5/theta-RCE_EXP_S50LN-AC-2012-01-01-000000-g1.h5' '/theta' }
        { 'HDF5/RCE_EXP_S50LN/HDF5/cloud-RCE_EXP_S50LN-AC-2012-01-01-000000-g1.h5' '/cloud' }
        { 'HDF5/RCE_EXP_S50LN/HDF5/vint_vapor-RCE_EXP_S50LN-AC-2012-01-01-000000-g1.h5' '/vertint_vapor' }
        'RCE_EXP_S50LN'
      }

      {
        { 'HDF5/RCE_EXP_S50LN_470/HDF5/tempc-RCE_EXP_S50LN_470-AC-2012-01-01-000000-g1.h5' '/tempc' }
        { 'HDF5/RCE_EXP_S50LN_470/HDF5/theta-RCE_EXP_S50LN_470-AC-2012-01-01-000000-g1.h5' '/theta' }
        { 'HDF5/RCE_EXP_S50LN_470/HDF5/cloud-RCE_EXP_S50LN_470-AC-2012-01-01-000000-g1.h5' '/cloud' }
        { 'HDF5/RCE_EXP_S50LN_470/HDF5/vint_vapor-RCE_EXP_S50LN_470-AC-2012-01-01-000000-g1.h5' '/vertint_vapor' }
        'RCE_EXP_S50LN_470'
      }

      {
        { 'HDF5/RCE_EXP_S70MY/HDF5/tempc-RCE_EXP_S70MY-AC-2012-01-01-000000-g1.h5' '/tempc' }
        { 'HDF5/RCE_EXP_S70MY/HDF5/theta-RCE_EXP_S70MY-AC-2012-01-01-000000-g1.h5' '/theta' }
        { 'HDF5/RCE_EXP_S70MY/HDF5/cloud-RCE_EXP_S70MY-AC-2012-01-01-000000-g1.h5' '/cloud' }
        { 'HDF5/RCE_EXP_S70MY/HDF5/vint_vapor-RCE_EXP_S70MY-AC-2012-01-01-000000-g1.h5' '/vertint_vapor' }
        'RCE_EXP_S70MY'
      }

%      {
%        { 'HDF5/RCE_BASE/HDF5/tempc-a-AC-2012-01-01-000000-g1.h5' '/tempc' }
%        { 'HDF5/RCE_BASE/HDF5/theta-a-AC-2012-01-01-000000-g1.h5' '/theta' }
%        { 'HDF5/RCE_BASE/HDF5/cloud-a-AC-2012-01-01-000000-g1.h5' '/cloud' }
%        { 'HDF5/RCE_BASE/HDF5/vint_vapor-a-AC-2012-01-01-000000-g1.h5' '/vertint_vapor' }
%        'RCE_BASE'
%      }

    };
  Nset = length(VarSets);

  PWthreshold = 40;


  fprintf('***************************************************************\n');
  fprintf('Generating vertical profile data:\n');

  for iset = 1:Nset
    TempFile   = VarSets{iset}{1}{1};
    TempVar    = VarSets{iset}{1}{2};

    ThetaFile   = VarSets{iset}{2}{1};
    ThetaVar    = VarSets{iset}{2}{2};

    CloudFile  = VarSets{iset}{3}{1};
    CloudVar   = VarSets{iset}{3}{2};

    PwFile  = VarSets{iset}{4}{1};
    PwVar   = VarSets{iset}{4}{2};

    Case       = VarSets{iset}{5};

    OutFile = sprintf('%s/SampleProfiles_%s.h5', Ddir, Case);

    fprintf('  Case: %s\n', Case);
    fprintf('\n');
    fprintf('    Precipitable water threshold: %f\n', PWthreshold);
    fprintf('\n');

    % Set up for extracting using nctools
    TEMP_DS = ncgeodataset(TempFile);
    THETA_DS = ncgeodataset(ThetaFile);
    CLOUD_DS = ncgeodataset(CloudFile);
    PW_DS = ncgeodataset(PwFile);

    TEMP_VAR = TEMP_DS.geovariable(TempVar);
    THETA_VAR = THETA_DS.geovariable(ThetaVar);
    CLOUD_VAR = CLOUD_DS.geovariable(CloudVar);
    PW_VAR = PW_DS.geovariable(PwVar);

    X_VAR = TEMP_DS.geovariable('x_coords');
    Y_VAR = TEMP_DS.geovariable('y_coords');
    Z_VAR = TEMP_DS.geovariable('z_coords');
    T_VAR = TEMP_DS.geovariable('t_coords');

    X = X_VAR.data(:);
    Y = Y_VAR.data(:);
    Z = Z_VAR.data(:);
    T = T_VAR.data(:);

    Nx = length(X);
    Ny = length(Y);
    Nz = length(Z);
    Nt = length(T);

    % Vars are organized as: (t,z,y,x)
    %   PW is (t,y,x)
    % Read in (y,x) data per time step and level
    % Then split (y,x) into two sets using PWthreshold
    % Form and average and record in output.
    fprintf('    Reading: %s (%s)\n', TempFile, TempVar);
    fprintf('    Reading: %s (%s)\n', ThetaFile, ThetaVar);
    fprintf('    Reading: %s (%s)\n', CloudFile, CloudVar);
    fprintf('    Reading: %s (%s)\n', PwFile, PwVar);
    fprintf('\n');

    % LAVG - average of low PW columns
    % HAVG - average of high PW columns
    TEMP_LAVG  = zeros([ Nz Nt ]);
    THETA_LAVG = zeros([ Nz Nt ]);
    CLOUD_LAVG = zeros([ Nz Nt ]);
    TEMP_HAVG  = zeros([ Nz Nt ]);
    THETA_HAVG = zeros([ Nz Nt ]);
    CLOUD_HAVG = zeros([ Nz Nt ]);

    % work one time step at a time and one level at a time
    % so that memory allocation does not get too big
    for it = 1:Nt
      % make PW a 3D with the (y,x) layers repeated on all z levels
      PW = squeeze(PW_VAR.data(it,:,:));

      for iz = 1:Nz
        TEMP  = squeeze(TEMP_VAR.data(it,iz,:,:));
        THETA = squeeze(THETA_VAR.data(it,iz,:,:));
        CLOUD = squeeze(CLOUD_VAR.data(it,iz,:,:));

        TEMP_LAVG(iz,it) = mean(TEMP(PW <= PWthreshold));
        TEMP_HAVG(iz,it) = mean(TEMP(PW >  PWthreshold));

        THETA_LAVG(iz,it) = mean(THETA(PW <= PWthreshold));
        THETA_HAVG(iz,it) = mean(THETA(PW >  PWthreshold));

        CLOUD_LAVG(iz,it) = mean(CLOUD(PW <= PWthreshold));
        CLOUD_HAVG(iz,it) = mean(CLOUD(PW >  PWthreshold));
      end

      if (mod(it,10) == 0)
        fprintf('    Working: time step = %d\n', it);
      end
    end
    fprintf('\n');

    % Write out the extracted profiles
    fprintf('    Writing: %s\n', OutFile);
    fprintf('\n');
    if (exist(OutFile, 'file') == 2)
      delete(OutFile);
    end

    OutX = 1;
    OutY = 1;
    OutZ = Z;
    OutT = T;

    OutNx = length(OutX);
    OutNy = length(OutY);
    OutNz = length(OutZ);
    OutNt = length(OutT);

    % Reshape the data to (x,y,z,t)
    % Vars are (z,t) at this point.
    TEMP_LAVG  = reshape(TEMP_LAVG,  [ OutNx OutNy OutNz OutNt ]);
    THETA_LAVG = reshape(THETA_LAVG, [ OutNx OutNy OutNz OutNt ]);
    CLOUD_LAVG = reshape(CLOUD_LAVG, [ OutNx OutNy OutNz OutNt ]);

    TEMP_HAVG  = reshape(TEMP_HAVG,  [ OutNx OutNy OutNz OutNt ]);
    THETA_HAVG = reshape(THETA_HAVG, [ OutNx OutNy OutNz OutNt ]);
    CLOUD_HAVG = reshape(CLOUD_HAVG, [ OutNx OutNy OutNz OutNt ]);

    TempLavgVname = sprintf('%s_lavg', TempVar);
    ThetaLavgVname = sprintf('%s_lavg', ThetaVar);
    CloudLavgVname = sprintf('%s_lavg', CloudVar);

    TempHavgVname = sprintf('%s_havg', TempVar);
    ThetaHavgVname = sprintf('%s_havg', ThetaVar);
    CloudHavgVname = sprintf('%s_havg', CloudVar);

    h5create(OutFile, TempLavgVname, size(TEMP_LAVG));
    h5create(OutFile, ThetaLavgVname, size(THETA_LAVG));
    h5create(OutFile, CloudLavgVname, size(CLOUD_LAVG));

    h5create(OutFile, TempHavgVname, size(TEMP_HAVG));
    h5create(OutFile, ThetaHavgVname, size(THETA_HAVG));
    h5create(OutFile, CloudHavgVname, size(CLOUD_HAVG));

    h5write(OutFile, TempLavgVname, TEMP_LAVG);
    h5write(OutFile, ThetaLavgVname, THETA_LAVG);
    h5write(OutFile, CloudLavgVname, CLOUD_LAVG);

    h5write(OutFile, TempHavgVname, TEMP_HAVG);
    h5write(OutFile, ThetaHavgVname, THETA_HAVG);
    h5write(OutFile, CloudHavgVname, CLOUD_HAVG);

    Xname = '/x_coords';
    Yname = '/y_coords';
    Zname = '/z_coords';
    Tname = '/t_coords';

    CreateDimensionsXyzt(OutFile, OutX, OutY, OutZ, OutT, Xname, Yname, Zname, Tname);
    AttachDimensionsXyzt(OutFile, TempLavgVname, Xname, Yname, Zname, Tname);
    AttachDimensionsXyzt(OutFile, ThetaLavgVname, Xname, Yname, Zname, Tname);
    AttachDimensionsXyzt(OutFile, CloudLavgVname, Xname, Yname, Zname, Tname);
    AttachDimensionsXyzt(OutFile, TempHavgVname, Xname, Yname, Zname, Tname);
    AttachDimensionsXyzt(OutFile, ThetaHavgVname, Xname, Yname, Zname, Tname);
    AttachDimensionsXyzt(OutFile, CloudHavgVname, Xname, Yname, Zname, Tname);
    NotateDimensionsXyzt(OutFile, Xname, Yname, Zname, Tname);
  end
end
