function [ ] = ExtractSampleProfiles()
% ExtractSampleProfiles extract vertical profile samples

  Ddir = 'DIAGS';

  % make sure output directory exists
  if (exist(Ddir, 'dir') ~= 7)
    mkdir(Ddir);
  end

  % input_file input_dataset case_name outfile_prefix output_dataset
  VarSets = {
%      {
%        { 'HDF5/RCE_EXP_S70LY/HDF5/tempc-RCE_EXP_S70LY-AC-2012-01-01-000000-g1.h5' '/tempc' }
%        { 'HDF5/RCE_EXP_S70LY/HDF5/theta-RCE_EXP_S70LY-AC-2012-01-01-000000-g1.h5' '/theta' }
%        { 'HDF5/RCE_EXP_S70LY/HDF5/cloud-RCE_EXP_S70LY-AC-2012-01-01-000000-g1.h5' '/cloud' }
%        { 'HDF5/RCE_EXP_S70LY/HDF5/vapor-RCE_EXP_S70LY-AC-2012-01-01-000000-g1.h5' '/vapor' }
%        'RCE_EXP_S70LY_cloudy'
%        182
%         13
%      }
%
%      {
%        { 'HDF5/RCE_EXP_S70LY/HDF5/tempc-RCE_EXP_S70LY-AC-2012-01-01-000000-g1.h5' '/tempc' }
%        { 'HDF5/RCE_EXP_S70LY/HDF5/theta-RCE_EXP_S70LY-AC-2012-01-01-000000-g1.h5' '/theta' }
%        { 'HDF5/RCE_EXP_S70LY/HDF5/cloud-RCE_EXP_S70LY-AC-2012-01-01-000000-g1.h5' '/cloud' }
%        { 'HDF5/RCE_EXP_S70LY/HDF5/vapor-RCE_EXP_S70LY-AC-2012-01-01-000000-g1.h5' '/vapor' }
%        'RCE_EXP_S70LY_clear'
%        1500
%         100
%      }
%
%      {
%        { 'HDF5/RCE_EXP_S70LN/HDF5/tempc-RCE_EXP_S70LN-AC-2012-01-01-000000-g1.h5' '/tempc' }
%        { 'HDF5/RCE_EXP_S70LN/HDF5/theta-RCE_EXP_S70LN-AC-2012-01-01-000000-g1.h5' '/theta' }
%        { 'HDF5/RCE_EXP_S70LN/HDF5/cloud-RCE_EXP_S70LN-AC-2012-01-01-000000-g1.h5' '/cloud' }
%        { 'HDF5/RCE_EXP_S70LN/HDF5/vapor-RCE_EXP_S70LN-AC-2012-01-01-000000-g1.h5' '/vapor' }
%        'RCE_EXP_S70LN_cloudy'
%        1303
%          46
%      }
%
%      {
%        { 'HDF5/RCE_EXP_S70LN/HDF5/tempc-RCE_EXP_S70LN-AC-2012-01-01-000000-g1.h5' '/tempc' }
%        { 'HDF5/RCE_EXP_S70LN/HDF5/theta-RCE_EXP_S70LN-AC-2012-01-01-000000-g1.h5' '/theta' }
%        { 'HDF5/RCE_EXP_S70LN/HDF5/cloud-RCE_EXP_S70LN-AC-2012-01-01-000000-g1.h5' '/cloud' }
%        { 'HDF5/RCE_EXP_S70LN/HDF5/vapor-RCE_EXP_S70LN-AC-2012-01-01-000000-g1.h5' '/vapor' }
%        'RCE_EXP_S70LN_clear'
%        936
%         43
%      }

      {
        { 'HDF5/RCE_EXP_S50LN/HDF5/tempc-RCE_EXP_S50LN-AC-2012-01-01-000000-g1.h5' '/tempc' }
        { 'HDF5/RCE_EXP_S50LN/HDF5/theta-RCE_EXP_S50LN-AC-2012-01-01-000000-g1.h5' '/theta' }
        { 'HDF5/RCE_EXP_S50LN/HDF5/cloud-RCE_EXP_S50LN-AC-2012-01-01-000000-g1.h5' '/cloud' }
        { 'HDF5/RCE_EXP_S50LN/HDF5/vapor-RCE_EXP_S50LN-AC-2012-01-01-000000-g1.h5' '/vapor' }
        'RCE_EXP_S50LN_cloudy'
        1986
          25
      }

      {
        { 'HDF5/RCE_EXP_S50LN/HDF5/tempc-RCE_EXP_S50LN-AC-2012-01-01-000000-g1.h5' '/tempc' }
        { 'HDF5/RCE_EXP_S50LN/HDF5/theta-RCE_EXP_S50LN-AC-2012-01-01-000000-g1.h5' '/theta' }
        { 'HDF5/RCE_EXP_S50LN/HDF5/cloud-RCE_EXP_S50LN-AC-2012-01-01-000000-g1.h5' '/cloud' }
        { 'HDF5/RCE_EXP_S50LN/HDF5/vapor-RCE_EXP_S50LN-AC-2012-01-01-000000-g1.h5' '/vapor' }
        'RCE_EXP_S50LN_clear'
        1500
         100
      }

    };
  Nset = length(VarSets);


  fprintf('***************************************************************\n');
  fprintf('Generating vertical profile data:\n');

  for iset = 1:Nset
    TempFile   = VarSets{iset}{1}{1};
    TempVar    = VarSets{iset}{1}{2};

    ThetaFile   = VarSets{iset}{2}{1};
    ThetaVar    = VarSets{iset}{2}{2};

    CloudFile  = VarSets{iset}{3}{1};
    CloudVar   = VarSets{iset}{3}{2};

    VaporFile  = VarSets{iset}{4}{1};
    VaporVar   = VarSets{iset}{4}{2};

    Case       = VarSets{iset}{5};
    Isamp      = VarSets{iset}{6};
    Jsamp      = VarSets{iset}{7};

    OutFile = sprintf('%s/SampleProfiles_%s.h5', Ddir, Case);

    fprintf('  Case: %s\n', Case);
    fprintf('\n');
    fprintf('    Extracting from (i,j) location: (%d,%d)\n', Isamp, Jsamp);
    fprintf('\n');

    % Set up for extracting using nctools
    TEMP_DS = ncgeodataset(TempFile);
    THETA_DS = ncgeodataset(ThetaFile);
    CLOUD_DS = ncgeodataset(CloudFile);
    VAPOR_DS = ncgeodataset(VaporFile);

    TEMP_VAR = TEMP_DS.geovariable(TempVar);
    THETA_VAR = THETA_DS.geovariable(ThetaVar);
    CLOUD_VAR = CLOUD_DS.geovariable(CloudVar);
    VAPOR_VAR = VAPOR_DS.geovariable(VaporVar);

    Z_VAR = TEMP_DS.geovariable('z_coords');
    T_VAR = TEMP_DS.geovariable('t_coords');

    % Vars are organized as: (t,z,y,x)
    % Extract a time series of the vertical profile located at (Isamp,Jsamp)
    % Take the specified column, plus the eight surrounding columns and form
    % an average profile from these.
    fprintf('    Reading: %s (%s)\n', TempFile, TempVar);
    TEMP = TEMP_VAR.data(:,:,Jsamp-1:Jsamp+1,Isamp-1:Isamp+1);

    fprintf('    Reading: %s (%s)\n', ThetaFile, ThetaVar);
    THETA = THETA_VAR.data(:,:,Jsamp-1:Jsamp+1,Isamp-1:Isamp+1);

    fprintf('    Reading: %s (%s)\n', CloudFile, CloudVar);
    CLOUD = CLOUD_VAR.data(:,:,Jsamp-1:Jsamp+1,Isamp-1:Isamp+1);
    
    fprintf('    Reading: %s (%s)\n', VaporFile, VaporVar);
    VAPOR = VAPOR_VAR.data(:,:,Jsamp-1:Jsamp+1,Isamp-1:Isamp+1);

    fprintf('\n');

    X = 1;
    Y = 1;
    Z = Z_VAR.data(:);
    T = T_VAR.data(:);

    Nx = 1;
    Ny = 1;
    Nz = length(Z);
    Nt = length(T);

    % Form the averages. Vars are: (t,z,y,x)
    TEMP_AVG = squeeze(mean(mean(TEMP,4),3));
    THETA_AVG = squeeze(mean(mean(THETA,4),3));
    CLOUD_AVG = squeeze(mean(mean(CLOUD,4),3));
    VAPOR_AVG = squeeze(mean(mean(VAPOR,4),3));

    % Sample specified column
    TEMP_SAMP = squeeze(TEMP(:,:,2,2));
    THETA_SAMP = squeeze(THETA(:,:,2,2));
    CLOUD_SAMP = squeeze(CLOUD(:,:,2,2));
    VAPOR_SAMP = squeeze(VAPOR(:,:,2,2));
    
    % Averaged profiles are organized as: (t,z)
    % Rearrange as (x,y,z,t) for subsequent use of ReadXyzt routine
    TEMP_AVG = permute(TEMP_AVG, [ 2 1 ]);
    TEMP_AVG = reshape(TEMP_AVG, [ Nx Ny Nz Nt ]);

    THETA_AVG = permute(THETA_AVG, [ 2 1 ]);
    THETA_AVG = reshape(THETA_AVG, [ Nx Ny Nz Nt ]);

    CLOUD_AVG = permute(CLOUD_AVG, [ 2 1 ]);
    CLOUD_AVG = reshape(CLOUD_AVG, [ Nx Ny Nz Nt ]);

    VAPOR_AVG = permute(VAPOR_AVG, [ 2 1 ]);
    VAPOR_AVG = reshape(VAPOR_AVG, [ Nx Ny Nz Nt ]);

    TEMP_SAMP = permute(TEMP_SAMP, [ 2 1 ]);
    TEMP_SAMP = reshape(TEMP_SAMP, [ Nx Ny Nz Nt ]);

    THETA_SAMP = permute(THETA_SAMP, [ 2 1 ]);
    THETA_SAMP = reshape(THETA_SAMP, [ Nx Ny Nz Nt ]);

    CLOUD_SAMP = permute(CLOUD_SAMP, [ 2 1 ]);
    CLOUD_SAMP = reshape(CLOUD_SAMP, [ Nx Ny Nz Nt ]);

    VAPOR_SAMP = permute(VAPOR_SAMP, [ 2 1 ]);
    VAPOR_SAMP = reshape(VAPOR_SAMP, [ Nx Ny Nz Nt ]);

    % Write out the extracted profiles
    fprintf('    Writing: %s\n', OutFile);
    fprintf('\n');
    if (exist(OutFile, 'file') == 2)
      delete(OutFile);
    end

    TempAvgVname = sprintf('%s_avg', TempVar);
    ThetaAvgVname = sprintf('%s_avg', ThetaVar);
    CloudAvgVname = sprintf('%s_avg', CloudVar);
    VaporAvgVname = sprintf('%s_avg', VaporVar);

    TempSampVname = sprintf('%s_samp', TempVar);
    ThetaSampVname = sprintf('%s_samp', ThetaVar);
    CloudSampVname = sprintf('%s_samp', CloudVar);
    VaporSampVname = sprintf('%s_samp', VaporVar);

    h5create(OutFile, TempAvgVname, size(TEMP_AVG));
    h5create(OutFile, ThetaAvgVname, size(THETA_AVG));
    h5create(OutFile, CloudAvgVname, size(CLOUD_AVG));
    h5create(OutFile, VaporAvgVname, size(VAPOR_AVG));

    h5create(OutFile, TempSampVname, size(TEMP_SAMP));
    h5create(OutFile, ThetaSampVname, size(THETA_SAMP));
    h5create(OutFile, CloudSampVname, size(CLOUD_SAMP));
    h5create(OutFile, VaporSampVname, size(VAPOR_SAMP));

    h5write(OutFile, TempAvgVname, TEMP_AVG);
    h5write(OutFile, ThetaAvgVname, THETA_AVG);
    h5write(OutFile, CloudAvgVname, CLOUD_AVG);
    h5write(OutFile, VaporAvgVname, VAPOR_AVG);

    h5write(OutFile, TempSampVname, TEMP_SAMP);
    h5write(OutFile, ThetaSampVname, THETA_SAMP);
    h5write(OutFile, CloudSampVname, CLOUD_SAMP);
    h5write(OutFile, VaporSampVname, VAPOR_SAMP);

    Xname = '/x_coords';
    Yname = '/y_coords';
    Zname = '/z_coords';
    Tname = '/t_coords';

    CreateDimensionsXyzt(OutFile, X, Y, Z, T, Xname, Yname, Zname, Tname);
    AttachDimensionsXyzt(OutFile, TempAvgVname, Xname, Yname, Zname, Tname);
    AttachDimensionsXyzt(OutFile, ThetaAvgVname, Xname, Yname, Zname, Tname);
    AttachDimensionsXyzt(OutFile, CloudAvgVname, Xname, Yname, Zname, Tname);
    AttachDimensionsXyzt(OutFile, VaporAvgVname, Xname, Yname, Zname, Tname);
    AttachDimensionsXyzt(OutFile, TempSampVname, Xname, Yname, Zname, Tname);
    AttachDimensionsXyzt(OutFile, ThetaSampVname, Xname, Yname, Zname, Tname);
    AttachDimensionsXyzt(OutFile, CloudSampVname, Xname, Yname, Zname, Tname);
    AttachDimensionsXyzt(OutFile, VaporSampVname, Xname, Yname, Zname, Tname);
    NotateDimensionsXyzt(OutFile, Xname, Yname, Zname, Tname);
  end
end
