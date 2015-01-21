function [ ] = ExtractSampleProfiles(ConfigFile)
% ExtractSampleProfiles extract vertical profile samples

  % Read the config file to get the structure of how the data is laid out in
  % the file system.
  [ Config ] = ReadConfig(ConfigFile);

  Ddir = Config.DiagDir;

  % make sure output directory exists
  if (exist(Ddir, 'dir') ~= 7)
    mkdir(Ddir);
  end

  % input_file input_dataset case_name outfile_prefix output_dataset
  VarSets = {
      {
        { 'HDF5/RCE50_RECT/tempc-a-AS-2012-01-01-000000-g1.h5' 'tempc' }
        { 'HDF5/RCE50_RECT/cloud-a-AS-2012-01-01-000000-g1.h5' 'cloud' }
        { 'HDF5/RCE50_RECT/vapor-a-AS-2012-01-01-000000-g1.h5' 'vapor' }
        'RCE50_RECT_800'
        800
        100
      }

      {
        { 'HDF5/RCE50_RECT/tempc-a-AS-2012-01-01-000000-g1.h5' 'tempc' }
        { 'HDF5/RCE50_RECT/cloud-a-AS-2012-01-01-000000-g1.h5' 'cloud' }
        { 'HDF5/RCE50_RECT/vapor-a-AS-2012-01-01-000000-g1.h5' 'vapor' }
        'RCE50_RECT_MID'
        1500
        100
      }

%      {
%        { 'HDF5/RCE50_OLD_RI/tempc-a-AS-2012-01-01-000000-g1.h5' 'tempc' }
%        { 'HDF5/RCE50_OLD_RI/cloud-a-AS-2012-01-01-000000-g1.h5' 'cloud' }
%        { 'HDF5/RCE50_OLD_RI/vapor-a-AS-2012-01-01-000000-g1.h5' 'vapor' }
%        'RCE50_OLD_RI_MID'
%        1500
%        100
%      }
%
%      {
%        { 'HDF5/RCE70_OLD_RI/tempc-a-AS-2012-01-01-000000-g1.h5' 'tempc' }
%        { 'HDF5/RCE70_OLD_RI/cloud-a-AS-2012-01-01-000000-g1.h5' 'cloud' }
%        { 'HDF5/RCE70_OLD_RI/vapor-a-AS-2012-01-01-000000-g1.h5' 'vapor' }
%        'RCE70_OLD_RI_MID'
%        1500
%        100
%      }
    };
  Nset = length(VarSets);


  fprintf('***************************************************************\n');
  fprintf('Generating vertical profile data:\n');

  for iset = 1:Nset
    TempFile   = VarSets{iset}{1}{1};
    TempVar    = VarSets{iset}{1}{2};

%    CloudFile  = VarSets{iset}{2}{1};
%    CloudVar   = VarSets{iset}{2}{2};

%    VaporFile  = VarSets{iset}{3}{1};
%    VaporVar   = VarSets{iset}{3}{2};

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
%    CLOUD_DS = ncgeodataset(CloudFile);
%    VAPOR_DS = ncgeodataset(VaporFile);

    TEMP_VAR = TEMP_DS.geovariable(TempVar);
%    CLOUD_VAR = CLOUD_DS.geovariable(CloudVar);
%    VAPOR_VAR = VAPOR_DS.geovariable(VaporVar);

    Z_VAR = TEMP_DS.geovariable('z_coords');
    T_VAR = TEMP_DS.geovariable('t_coords');

    % Vars are organized as: (t,z,y,x)
    % Extract a time series of the vertical profile located at (Isamp,Jsamp)
    fprintf('    Reading: %s (%s)\n', TempFile, TempVar);
    TEMP = TEMP_VAR.data(:,:,Jsamp,Isamp);

%    fprintf('    Reading: %s (%s)\n', CloudFile, CloudVar);
%    CLOUD = CLOUD_VAR.data(:,:,Jsamp,Isamp);
    
%    fprintf('    Reading: %s (%s)\n', VaporFile, VaporVar);
%    VAPOR = VAPOR_VAR.data(:,:,Jsamp,Isamp);

    Z = Z_VAR.data(:);
    T = T_VAR.data(:);

    fprintf('\n');
    
    % Extracted profiles are organized as: (t,z)
    % Rearrange as (x,y,z,t) for subsequent use of ReadXyzt routine
    Nx = 1;
    Ny = 1;
    [ Nt, Nz ] = size(TEMP);

    TEMP = permute(TEMP, [ 2 1 ]);
    TEMP = reshape(TEMP, [ Nx Ny Nz Nt ]);

%    CLOUD = permute(CLOUD, [ 2 1 ]);
%    CLOUD = reshape(CLOUD, [ Nx Ny Nz Nt ]);

%    VAPOR = permute(VAPOR, [ 2 1 ]);
%    VAPOR = reshape(VAPOR, [ Nx Ny Nz Nt ]);

    % Write out the extracted profiles
    fprintf('    Writing: %s\n', OutFile);
    fprintf('\n');

    Xdummy = 1;
    Ydummy = 1;

    hdf5write(OutFile, TempVar, TEMP); 
%    hdf5write(OutFile, CloudVar, CLOUD, 'WriteMode', 'append'); 
%    hdf5write(OutFile, VaporVar, VAPOR, 'WriteMode', 'append'); 

    hdf5write(OutFile, '/x_coords', Xdummy, 'WriteMode', 'append');
    hdf5write(OutFile, '/y_coords', Ydummy, 'WriteMode', 'append');
    hdf5write(OutFile, '/z_coords', Z,      'WriteMode', 'append');
    hdf5write(OutFile, '/t_coords', T,      'WriteMode', 'append');
  end
end
