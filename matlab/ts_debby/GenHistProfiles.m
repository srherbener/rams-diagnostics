function [ ] = GenHistProfiles(ConfigFile)

  [ Config ] = ReadConfig(ConfigFile);

  Adir = Config.AzavgDir;
  Ddir = Config.DiagDir;

  % make sure output directory exists
  if (exist(Ddir, 'dir') ~= 7)
    mkdir(Ddir);
  end

  % Description of profiles
  %   <file_prefix> <var_name> <reduction_method> <arg_for_reduction_method> <out_var_name> <select_op> <select_val>
  ProfList = {
    { 'hist_ccn_conc' 'ccn_conc' 'wtmean' 0.5 'prof_ccn_conc'      'ge' 20   }
    { 'hist_d1_conc'  'd1_conc'  'wtmean' 0.5 'prof_d1_conc'       'ge' 20   }
    { 'hist_d2_conc'  'd2_conc'  'wtmean' 0.5 'prof_d2_conc'       'ge' 20   }

    { 'hist_lh_vapt'  'lh_vapt'  'wtmean' 0.5 'prof_lat_heat_vapt' 'ge'  10  }
    { 'hist_lh_vapt'  'lh_vapt'  'wtmean' 0.5 'prof_lat_cool_vapt' 'le' -10  }

    { 'hist_lh_frzt'  'lh_frzt'  'wtmean' 0.5 'prof_lat_heat_frzt' 'ge'  1   }
    { 'hist_lh_frzt'  'lh_frzt'  'wtmean' 0.5 'prof_lat_cool_frzt' 'le' -1   }

    { 'hist_w'        'w'        'wtmean' 0.5 'prof_w_up'          'ge'  0.1 }
    { 'hist_w'        'w'        'wtmean' 0.5 'prof_w_down'        'le' -0.1 }

    { 'hist_theta_e'  'theta_e'  'wtmean' 0.5 'prof_theta_e'       'ge'  0   }
    };

    Nprof = length(ProfList);

  for icase = 1:length(Config.Cases)
    Case = Config.Cases(icase).Cname;

    fprintf('*****************************************************************\n');
    fprintf('Generating profiles:\n');
    fprintf('  Case: %s\n', Case);
    fprintf('\n');

    % Put all profiles into one file per case
    % If the file exists, remove it so that the HDF5 commands
    % can effectively re-create datasets.
    OutFile = sprintf('%s/hist_profs_%s.h5', Ddir, Case);
    if (exist(OutFile, 'file') == 2)
      delete(OutFile);
    end

    for iprof = 1:Nprof
      Fprefix   = ProfList{iprof}{1};
      Vname     = ProfList{iprof}{2};
      Rmethod   = ProfList{iprof}{3};
      Ptile     = ProfList{iprof}{4};
      OutVname  = ProfList{iprof}{5};
      SelectOp  = ProfList{iprof}{6};
      SelectVal = ProfList{iprof}{7};

      % add on leading '/' for HDF5 routines
      Vname    = sprintf('/%s', Vname);
      OutVname = sprintf('/%s', OutVname);

      InFile = sprintf('%s/%s_%s.h5', Adir, Fprefix, Case);

      fprintf('  Reading: %s (%s)\n', InFile, Vname);
      fprintf('    Reduction method: %s (%.2f)\n', Rmethod, Ptile);
      fprintf('    Selection: %s %.2f\n', SelectOp, SelectVal);

      % Read in data which will be 4D -> (x,y,z,t)
      %
      %     x --> radial bands
      %     y --> histogram bins
      %     z --> height
      %     t --> time
      %
      HDATA = squeeze(h5read(InFile, Vname));
      BINS  = squeeze(h5read(InFile, '/y_coords'));

      % Assume same r,z,t values for all profiles
      if (iprof == 1)
        X    = squeeze(h5read(InFile, '/x_coords'));
        Y    = 1; % dummy dimension for output
        Z    = squeeze(h5read(InFile, '/z_coords'));
        T    = squeeze(h5read(InFile, '/t_coords'));

        Nx = length(X);
        Ny = length(Y);
        Nz = length(Z);
        Nt = length(T);
      end

      % Reduce the histograms level by level to create vertical profiles
      % PROF will be organized as: (r,z,t)
      %
      % Do selection on bin values (ge or le)
      %

      % select all bin values by default
      B1 = 1;         
      B2 = length(BINS);

      if (strcmp(SelectOp, 'ge'))
        B1 = find(BINS >= SelectVal, 1, 'first');
      end

      if (strcmp(SelectOp, 'le'))
        B2 = find(BINS <= SelectVal, 1, 'last');
      end

      PROF = squeeze(ReduceHists(HDATA(:,B1:B2,:,:), 2, BINS(B1:B2), Rmethod, Ptile));

      % Write out measurement
      fprintf('  Writing: %s (%s)\n', OutFile, OutVname)

      % Write out profile -> force to be (x,y,z,t) for dimension
      % attach code below.
      OutVar = reshape(PROF, [ Nx Ny Nz Nt ]);
      h5create(OutFile, OutVname, size(OutVar));
      h5write(OutFile, OutVname, OutVar);
    end

    % Create the dimensions
    fprintf('  Creating dimensions: %s\n', OutFile);
    Xname = '/x_coords';
    Yname = '/y_coords';
    Zname = '/z_coords';
    Tname = '/t_coords';

    CreateDimensionsXyzt(OutFile, X, Y, Z, T, Xname, Yname, Zname, Tname);

    % Attach dimensions to all variables
    for iprof = 1:Nprof
      Vname = ProfList{iprof}{5};   % use the output var name
      AttachDimensionsXyzt(OutFile, Vname, Xname, Yname, Zname, Tname);
    end

    % GRADS needs the following attributes on the dimension datasets in order
    % to recognize which dimensions go with which datasets. These attribute
    % names and values are following the COARDS convention.
    WriteStringAttribute(OutFile, Xname, 'axis', 'x');
    WriteStringAttribute(OutFile, Xname, 'units', 'degrees_east');

    WriteStringAttribute(OutFile, Yname, 'axis', 'y');
    WriteStringAttribute(OutFile, Yname, 'units', 'degrees_north');

    WriteStringAttribute(OutFile, Zname, 'axis', 'z');
    WriteStringAttribute(OutFile, Zname, 'units', 'meters');

    WriteStringAttribute(OutFile, Tname, 'axis', 't');
    WriteStringAttribute(OutFile, Tname, 'units', 'seconds since 2006-08-20 12:00:00 00:00');
    
    fprintf('\n');
  end
end 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CreateDimsensionsXyzt()
%
% This routine will create 4D dimensions in the HDF5 file.
%
%
function [] = CreateDimensionsXyzt(File, X, Y, Z, T, Xname, Yname, Zname, Tname)

  Nx = length(X);
  Ny = length(Y);
  Nz = length(Z);
  Nt = length(T);

  % write out the coordinate values
  h5create(File, Xname, Nx);
  h5write(File,  Xname, X);

  h5create(File, Yname, Ny);
  h5write(File,  Yname, Y);

  h5create(File, Zname, Nz);
  h5write(File,  Zname, Z);

  h5create(File, Tname, Nt);
  h5write(File,  Tname, T);

  % Use low level HDF5 routiens to mark these as dimensions
  file_id = H5F.open(File, 'H5F_ACC_RDWR', 'H5P_DEFAULT');

  x_id = H5D.open(file_id, Xname, 'H5P_DEFAULT');
  y_id = H5D.open(file_id, Yname, 'H5P_DEFAULT');
  z_id = H5D.open(file_id, Zname, 'H5P_DEFAULT');
  t_id = H5D.open(file_id, Tname, 'H5P_DEFAULT');

  % set x,y,z,t as dimensions
  H5DS.set_scale(x_id, 'x');
  H5DS.set_scale(y_id, 'y');
  H5DS.set_scale(z_id, 'z');
  H5DS.set_scale(t_id, 't');
  
  % close dimension variables, and file
  H5D.close(x_id);
  H5D.close(y_id);
  H5D.close(z_id);
  H5D.close(t_id);
  H5F.close(file_id);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WriteStringAttribute()
%
% This routine will create and attach a string attribute
% to the given file and dataset.
%
function [] = WriteStringAttribute(File, Dataset, AttrName, AttrString)

  file_id = H5F.open(File, 'H5F_ACC_RDWR', 'H5P_DEFAULT');
  dset_id = H5D.open(file_id, Dataset, 'H5P_DEFAULT'); 

  % Going to use C style strings which are null terminated, so append
  % the null character on the end of AttrString
  String = sprintf('%s%c', AttrString, char(0));
  Slen = length(String);

  atype_id = H5T.copy('H5T_C_S1');
  H5T.set_size(atype_id, Slen);

  % Use a scalar memory space so that a single string is created
  % (as opposed to an array of characters).
  mspc_id = H5S.create('H5S_SCALAR');

  attr_id = H5A.create(dset_id, AttrName, atype_id, mspc_id , 'H5P_DEFAULT');
  H5A.write(attr_id, atype_id, String);

  H5S.close(mspc_id);
  H5A.close(attr_id);
  H5T.close(atype_id);
  H5D.close(dset_id);
  H5F.close(file_id);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AttachDimensionsXyzt()
%
% This routine will attach the x,y,z,t dimensions to the
% given file and dataset. It is assumed that the dataset
% is organized as (x,y,z,t).
%
function [] = AttachDimensionsXyzt(File, Dataset, Xname, Yname, Zname, Tname);

  file_id = H5F.open(File, 'H5F_ACC_RDWR', 'H5P_DEFAULT');
  dset_id = H5D.open(file_id, Dataset, 'H5P_DEFAULT'); 

  x_id = H5D.open(file_id, Xname, 'H5P_DEFAULT');
  y_id = H5D.open(file_id, Yname, 'H5P_DEFAULT');
  z_id = H5D.open(file_id, Zname, 'H5P_DEFAULT');
  t_id = H5D.open(file_id, Tname, 'H5P_DEFAULT');

  % Dimensions go in reverse order since the data gets stored that
  % way. MATLAB is column-major, whereas C is used to write the HDF5
  % file making the file storage row-major.
  %
  H5DS.attach_scale(dset_id, x_id, 3);
  H5DS.attach_scale(dset_id, y_id, 2);
  H5DS.attach_scale(dset_id, z_id, 1);
  H5DS.attach_scale(dset_id, t_id, 0);

  H5D.close(x_id);
  H5D.close(y_id);
  H5D.close(z_id);
  H5D.close(t_id);
  H5D.close(dset_id);
  H5F.close(file_id);
end
