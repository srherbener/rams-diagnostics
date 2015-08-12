function [ ] = ExtractSampleProfile()
%ExtractSampleProfile extract a time series of sampled vertical profiles

  CaseList = {
        'TSD_SAL_DUST_TR'
        'TSD_NONSAL_DUST'
        };
       Nc = length(CaseList);
    
      %  This routine walk through the cases (CaseList) and samples (SampleList), read in
      %  the samples and select columns based on the arguments Xselect and Yselect.
      %  Then the columns are averaged into a single column and the time series of these
      %  averaged columms are output into the dataset OutVname inside OutFile.
      %
      %  Xselect and Yselect are strings specifying index ranges. Note that if the input dataset
      %  is 4D -> (x,y,z,t), then:
      %
      %    Xselect     Yselect          Selection
      %     '5'         '7'           Var(5,7,:,:)
      %    '1:3'      '4:2:10'        Var(1:3,4:2:10,:,:)
      %
    
      % InFile InVname Xselect Yselect OutFile OutVname
      SampleList = {
        { 'HDF5/<CASE>/HDF5/d1_mass-<CASE>-AS-2006-08-20-120000-g3.h5' '/dust1_mass' '133:135' '468:470' 'DIAGS/sp_dust1_<CASE>.h5' '/sp_dust1_mass' }
        { 'HDF5/<CASE>/HDF5/d2_mass-<CASE>-AS-2006-08-20-120000-g3.h5' '/dust2_mass' '133:135' '468:470' 'DIAGS/sp_dust2_<CASE>.h5' '/sp_dust2_mass' }
        { 'HDF5/<CASE>/HDF5/tr3d1-<CASE>-AS-2006-08-20-120000-g3.h5' '/tr3d1' '133:135' '468:470' 'DIAGS/sp_tracer1_<CASE>.h5' '/sp_tr3d1' }
        { 'HDF5/<CASE>/HDF5/tr3d2-<CASE>-AS-2006-08-20-120000-g3.h5' '/tr3d2' '133:135' '468:470' 'DIAGS/sp_tracer2_<CASE>.h5' '/sp_tr3d2' }
        };
      Ns = length(SampleList);

      for ic = 1:Nc
        Case = CaseList{ic};
        fprintf('*************************************************************************\n');
        fprintf('Extracting sample profile for case: %s\n', Case);
    
        for is = 1:Ns
          InFile   = regexprep(SampleList{is}{1}, '<CASE>', Case);
          InVname  = SampleList{is}{2};
          Xselect  = SampleList{is}{3};
          Yselect  = SampleList{is}{4};
          OutFile  = regexprep(SampleList{is}{5}, '<CASE>', Case);
          OutVname = SampleList{is}{6};
    
          fprintf('  Reading: %s (%s)\n', InFile, InVname);
          fprintf('    Xselect: %s\n', Xselect);
          fprintf('    Yselect: %s\n', Yselect);
          fprintf('\n');
         
          IN_DS  = ncgeodataset(InFile);
        
          IN_VAR = IN_DS.geovariable(InVname);
          X_VAR  = IN_DS.geovariable('/x_coords');
          Y_VAR  = IN_DS.geovariable('/y_coords');
          Z_VAR  = IN_DS.geovariable('/z_coords');
          T_VAR  = IN_DS.geovariable('/t_coords');
        
          S = IN_VAR.size;
          Ndims = length(S);
        
          % Assume dimensionality of IN_VAR is one of:
          %   Ndims == 1 -> (t)
          %   Ndims == 2 -> (t,z)
          %   Ndims == 3 -> (t,z,y)
          %   Ndims == 4 -> (t,z,y,x)
        
          Nt = S(1);
          Nz = 1;
          if (Ndims > 1)
            Nz = S(2);
          end
        
          % Walk through each time step. Read in the var, and if needed compute the mean
          % across the x and y ranges. The input dataset will be a time series of either
          % 1D, 2D or 3D fields. Assume the following cases depending on the number of
          % dimensions.
          %
          %  1D -> (t,z)       Ndims == 2
          %  2D -> (t,z,y)     Ndims == 3
          %  3D -> (t,z,y,x)   Ndims == 4
          %
          % This means:
          %
          %  Ndims == 1 --> read produces a single point
          %
          %  Ndims == 2 --> read produces 1D: (z)
          %
          %  Ndims == 3 --> select y on read
          %                 read produces 2D: (z,y),
          %                 take mean on dimension 2
          %
          %  Ndims == 4 --> select y and x on read
          %                 read produces 3D: (z,y,x),
          %                 take mean on dimensions 2 and 3
          %
        
          PROF_TS = squeeze(zeros([ Nz Nt ]));
          for it = 1:Nt
            % meassage so user can see that progress is being made
            if (mod(it,10) == 0)
              fprintf('    Working, time step = %d\n', it);
            end
        
            % set up the read and select command
            if (Ndims == 1)
              ReadSelect = sprintf('squeeze(IN_VAR.data(%d))', it);
            elseif (Ndims == 2)
              ReadSelect = sprintf('squeeze(IN_VAR.data(%d,:))', it);
            elseif (Ndims == 3)
              ReadSelect = sprintf('squeeze(IN_VAR.data(%d,:,%s))', it, Yselect);
            elseif (Ndims == 4)
              ReadSelect = sprintf('squeeze(IN_VAR.data(%d,:,%s,%s))', it, Yselect, Xselect);
            end
        
            % do the read, select
            VAR = eval(ReadSelect);
        
            % if necessary do use nanmean to reduce to a single profile
            for i = Ndims-1:-1:2
              VAR = squeeze(nanmean(VAR,i));
            end
        
            % Var has been reduced to either a single point or a single profile
            PROF_TS(:,it) = VAR;
          end
          fprintf('\n');
        
        
          % Write output so that grads can read in the file:
          %  4D -> (x,y,z,t), and x and y are dummy variables.
        
          % Set coordinates, X and Y are dummy values since these have been
          % eliminated by the averaging.
          X = 1; 
          Y = 1;
          Z = Z_VAR.data(:);
          T = T_VAR.data(:);
        
          Nx = 1;
          Ny = 1;
        
          % PROF_TS is (z,t)
          OutVar = reshape(PROF_TS, [ Nx Ny Nz Nt] );
        
          fprintf('  Writing: %s (%s)\n', OutFile, OutVname);
          if (exist(OutFile, 'file') == 2)
            delete(OutFile);
          end
        
          h5create(OutFile, OutVname, size(OutVar));
          h5write (OutFile, OutVname, OutVar);
        
          Xname = '/x_coords';
          Yname = '/y_coords';
          Zname = '/z_coords';
          Tname = '/t_coords';
          
          CreateDimensionsXyzt(OutFile, X, Y, Z, T, Xname, Yname, Zname, Tname);
          AttachDimensionsXyzt(OutFile, OutVname, Xname, Yname, Zname, Tname);
          NotateDimensionsXyzt(OutFile, Xname, Yname, Zname, Tname);

          fprintf('\n');
        end
      end
    

end
