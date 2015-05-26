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
        { 'HDF5/RCE_CNTL/HDF5/tempc-RCE_CNTL-AC-2012-01-01-000000-g1.h5' '/tempc' }
        { 'HDF5/RCE_CNTL/HDF5/cloud-RCE_CNTL-AC-2012-01-01-000000-g1.h5' '/cloud' }
        { 'HDF5/RCE_CNTL/HDF5/vapor-RCE_CNTL-AC-2012-01-01-000000-g1.h5' '/vapor' }
        'RCE_CNTL_NEW_DIFF_cloudy'
        182
         13
      }

      {
        { 'HDF5/RCE_CNTL/HDF5/tempc-RCE_CNTL-AC-2012-01-01-000000-g1.h5' '/tempc' }
        { 'HDF5/RCE_CNTL/HDF5/cloud-RCE_CNTL-AC-2012-01-01-000000-g1.h5' '/cloud' }
        { 'HDF5/RCE_CNTL/HDF5/vapor-RCE_CNTL-AC-2012-01-01-000000-g1.h5' '/vapor' }
        'RCE_CNTL_NEW_DIFF_clear'
        1500
         100
      }

      {
        { 'HDF5/RCE_CNTL/HDF5_NZ75/tempc-RCE_CNTL-AC-2012-01-01-000000-g1.h5' '/tempc' }
        { 'HDF5/RCE_CNTL/HDF5_NZ75/cloud-RCE_CNTL-AC-2012-01-01-000000-g1.h5' '/cloud' }
        { 'HDF5/RCE_CNTL/HDF5_NZ75/vapor-RCE_CNTL-AC-2012-01-01-000000-g1.h5' '/vapor' }
        'RCE_CNTL_OLD_DIFF_cloudy'
        1303
          46
      }

      {
        { 'HDF5/RCE_CNTL/HDF5_NZ75/tempc-RCE_CNTL-AC-2012-01-01-000000-g1.h5' '/tempc' }
        { 'HDF5/RCE_CNTL/HDF5_NZ75/cloud-RCE_CNTL-AC-2012-01-01-000000-g1.h5' '/cloud' }
        { 'HDF5/RCE_CNTL/HDF5_NZ75/vapor-RCE_CNTL-AC-2012-01-01-000000-g1.h5' '/vapor' }
        'RCE_CNTL_OLD_DIFF_clear'
        936
         43
      }

    };
  Nset = length(VarSets);


  fprintf('***************************************************************\n');
  fprintf('Generating vertical profile data:\n');

  for iset = 1:Nset
    TempFile   = VarSets{iset}{1}{1};
    TempVar    = VarSets{iset}{1}{2};

    CloudFile  = VarSets{iset}{2}{1};
    CloudVar   = VarSets{iset}{2}{2};

    VaporFile  = VarSets{iset}{3}{1};
    VaporVar   = VarSets{iset}{3}{2};

    Case       = VarSets{iset}{4};
    Isamp      = VarSets{iset}{5};
    Jsamp      = VarSets{iset}{6};

    OutFile = sprintf('%s/SampleProfiles_%s.h5', Ddir, Case);

    fprintf('  Case: %s\n', Case);
    fprintf('\n');
    fprintf('    Extracting from (i,j) location: (%d,%d)\n', Isamp, Jsamp);
    fprintf('\n');

    % Set up for extracting using nctools
    TEMP_DS = ncgeodataset(TempFile);
    CLOUD_DS = ncgeodataset(CloudFile);
    VAPOR_DS = ncgeodataset(VaporFile);

    TEMP_VAR = TEMP_DS.geovariable(TempVar);
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
    CLOUD_AVG = squeeze(mean(mean(CLOUD,4),3));
    VAPOR_AVG = squeeze(mean(mean(VAPOR,4),3));
    
    % Averaged profiles are organized as: (t,z)
    % Rearrange as (x,y,z,t) for subsequent use of ReadXyzt routine
    TEMP_AVG = permute(TEMP_AVG, [ 2 1 ]);
    TEMP_AVG = reshape(TEMP_AVG, [ Nx Ny Nz Nt ]);

    CLOUD_AVG = permute(CLOUD_AVG, [ 2 1 ]);
    CLOUD_AVG = reshape(CLOUD_AVG, [ Nx Ny Nz Nt ]);

    VAPOR_AVG = permute(VAPOR_AVG, [ 2 1 ]);
    VAPOR_AVG = reshape(VAPOR_AVG, [ Nx Ny Nz Nt ]);

    % Write out the extracted profiles
    fprintf('    Writing: %s\n', OutFile);
    fprintf('\n');
    if (exist(OutFile, 'file') == 2)
      delete(OutFile);
    end

    h5create(OutFile, TempVar, size(TEMP_AVG));
    h5create(OutFile, CloudVar, size(CLOUD_AVG));
    h5create(OutFile, VaporVar, size(VAPOR_AVG));

    h5write(OutFile, TempVar, TEMP_AVG);
    h5write(OutFile, CloudVar, CLOUD_AVG);
    h5write(OutFile, VaporVar, VAPOR_AVG);

    Xname = '/x_coords';
    Yname = '/y_coords';
    Zname = '/z_coords';
    Tname = '/t_coords';

    CreateDimensionsXyzt(OutFile, X, Y, Z, T, Xname, Yname, Zname, Tname);
    AttachDimensionsXyzt(OutFile, TempVar, Xname, Yname, Zname, Tname);
    AttachDimensionsXyzt(OutFile, CloudVar, Xname, Yname, Zname, Tname);
    AttachDimensionsXyzt(OutFile, VaporVar, Xname, Yname, Zname, Tname);
    NotateDimensionsXyzt(OutFile, Xname, Yname, Zname, Tname);
  end
end
