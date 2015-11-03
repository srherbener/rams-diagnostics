function [ ] = GenRceHistMeas()

  % make sure output directory exists
  Ddir = 'DIAGS';  % coordinate this with output file names below
  if (exist(Ddir, 'dir') ~= 7)
    mkdir(Ddir);
  end

  % list of simulation cases
  CaseList = {
    'RCE_S300_CONV'
    };
  Ncases = length(CaseList);

  % Temporal ranges
  TstartPreSal = 10;
  TendPreSal   = 30;

  TstartSal = 40;
  TendSal   = 60;

  % Input is a histogram based upon either 2D or 3D field
  %
  % Azimuthal averaged data (azavg output):
  %   2D - input will be of the form (r,b,t)
  %   3D - input will be of the form (r,b,z,t)
  %
  % Time series averaged data (tsavg output):
  %   2D - input will be of the form (b,t)
  %   3D - input will be of the form (b,z,t)
  %
  %   where: r - radius
  %          b - bins (histogram counts)
  %          z - height
  %          t - time

  % The input form and reduction is specified by a 4 letter string. The string follows
  % the convention of azavg and tsavg output always writing a 4 dimension variable of
  % the form (x,y,z,t). The reduction spec uses 4 letters that correspond to the
  % xyzt dimension order of the input. If a lower case string is used, then that means
  % to reduce that dimension. If an underscore is used, that means that this input
  % dimension (one of x,y,z,t) is not used. For example, to read in a 3D field from azavg, and reduce
  % to a height profile, use 'rbZt'. This means:
  %
  %   input dim     quatity     reduce
  %       x          radius       yes
  %       y          bins         yes
  %       z          height       no
  %       t          time         yes
  %
  % To produce a time series of height-radius cross sections from azavg 3D output,
  % use 'RbZT'. This means:
  %
  %   input dim     quatity     reduce
  %       x          radius       no
  %       y          bins         yes
  %       z          height       no
  %       t          time         no
  %
  % To produce a radial profile from a 2D azavg output, use 'Rb_t'. This means:
  %
  %
  %   input dim     quatity     reduce
  %       x          radius       no
  %       y          bins         yes
  %       z          unused       ---
  %       t          time         yes
  %
  %
  % The idea for reduction is to combine (sum) all of the histogram counts first, then
  % reduce the histogram counts once (as opposed to reducing the histogram counts first
  % followed by averaging to combine along the other dimensions).
  %
  % Another notion is to avoid combining counts along the z axis, i.e., preserve
  % levels. This means that to get an r_profile from a 3D field, a single level
  % must be selected. 
  %
  % The temporal, radial and height ranges above are used for reduction. Temporal and
  % radial ranges are used for summing histogram counts, and height is used for selecting
  % a single level.
  %
  % Range specs for summing counts:
  %
  %    Radial range
  %      core            RstartCore,  RendCore
  %      rband           RstartRband, RendRband
  %      env             RstartEnv,   RendEnv
  %
  %    Temporal range
  %       pre_sal        TstartPreSal, TendPreSal
  %       sal            TstartSal,    TendSal
  %
  % Range specs for selecting a level:
  %
  %    Height range
  %       sfc            Zsfc
  %       maxlev         Go through histogram reduction steps, then select level that
  %                      contains the maximum value AND lies between Zsfc and Ztop.
  %       minlev         Same as maxlev, except select the level with the minumum value.

  % Description of measurements
  %   <measure_name> <measure_list> <out_file>
  %
  %   where <measure_list> is one or more of:
  %     <in_file> <var_name> <reduction_method> <arg_for_reduction_method> <out_var_name> <reduce_spec> <r_range> <t_range> <z_range> <bin_select_op> <bin_select_val>
  %       ranges are defined above
  %       <reduce_spec> describes input dimensions (and order) with upper case letters on the dims that are to be reduced.
  %           <reduce_spec> = 4-character string that corresponds to 'xyzt' that describes which dimensions are radius, bins, height and time,
  %                           as well as sepcifiying which dimensions are to be reduced (lower case).
  %                For example, 'rb_T' means that the input has three dimensions (after squeezing), radius is x, bins are y, and t (z will disappaer
  %                after squeezing), making the input data (r,b,t). The lower case 'r' and 'b' say to reduce the r and b dimensions, leaving
  %                the output a vector in dimension t.
  %
  %                'RbZT' means radius is x, bins are y, height is z and time is t, reduce only the bins, leaving the output (r,z,t).

  MeasSets = {

    {
      'Cold Pools Tsavg'
      {
        % high precip regions
        { 'HDF5/TsAveragedData/hist_lhf_cool_hp_<CASE>.h5' '/hist_lhf_cool' 'wtmean'  0.0  '/lhf_cool_hp_ts'      'b_ZT'   ''      ''        ''       'le' 0 }
        { 'HDF5/TsAveragedData/hist_lhf_cool_hp_<CASE>.h5' '/hist_lhf_cool' 'wtmean'  0.0  '/lhf_cool_hp'         'b_Zt'   ''      ''        ''       'le' 0 }

        { 'HDF5/TsAveragedData/hist_lhf_heat_hp_<CASE>.h5' '/hist_lhf_heat' 'wtmean'  0.0  '/lhf_heat_hp_ts'      'b_ZT'   ''      ''        ''       'ge' 0 }
        { 'HDF5/TsAveragedData/hist_lhf_heat_hp_<CASE>.h5' '/hist_lhf_heat' 'wtmean'  0.0  '/lhf_heat_hp'         'b_Zt'   ''      ''        ''       'ge' 0 }

        { 'HDF5/TsAveragedData/hist_lhv_cool_hp_<CASE>.h5' '/hist_lhv_cool' 'wtmean'  0.0  '/lhv_cool_hp_ts'      'b_ZT'   ''      ''        ''       'le' 0 }
        { 'HDF5/TsAveragedData/hist_lhv_cool_hp_<CASE>.h5' '/hist_lhv_cool' 'wtmean'  0.0  '/lhv_cool_hp'         'b_Zt'   ''      ''        ''       'le' 0 }

        { 'HDF5/TsAveragedData/hist_lhv_heat_hp_<CASE>.h5' '/hist_lhv_heat' 'wtmean'  0.0  '/lhv_heat_hp_ts'      'b_ZT'   ''      ''        ''       'ge' 0 }
        { 'HDF5/TsAveragedData/hist_lhv_heat_hp_<CASE>.h5' '/hist_lhv_heat' 'wtmean'  0.0  '/lhv_heat_hp'         'b_Zt'   ''      ''        ''       'ge' 0 }

        % low precip regions
        { 'HDF5/TsAveragedData/hist_lhf_cool_lp_<CASE>.h5' '/hist_lhf_cool' 'wtmean'  0.0  '/lhf_cool_lp_ts'      'b_ZT'   ''      ''        ''       'le' 0 }
        { 'HDF5/TsAveragedData/hist_lhf_cool_lp_<CASE>.h5' '/hist_lhf_cool' 'wtmean'  0.0  '/lhf_cool_lp'         'b_Zt'   ''      ''        ''       'le' 0 }

        { 'HDF5/TsAveragedData/hist_lhf_heat_lp_<CASE>.h5' '/hist_lhf_heat' 'wtmean'  0.0  '/lhf_heat_lp_ts'      'b_ZT'   ''      ''        ''       'ge' 0 }
        { 'HDF5/TsAveragedData/hist_lhf_heat_lp_<CASE>.h5' '/hist_lhf_heat' 'wtmean'  0.0  '/lhf_heat_lp'         'b_Zt'   ''      ''        ''       'ge' 0 }

        { 'HDF5/TsAveragedData/hist_lhv_cool_lp_<CASE>.h5' '/hist_lhv_cool' 'wtmean'  0.0  '/lhv_cool_lp_ts'      'b_ZT'   ''      ''        ''       'le' 0 }
        { 'HDF5/TsAveragedData/hist_lhv_cool_lp_<CASE>.h5' '/hist_lhv_cool' 'wtmean'  0.0  '/lhv_cool_lp'         'b_Zt'   ''      ''        ''       'le' 0 }

        { 'HDF5/TsAveragedData/hist_lhv_heat_lp_<CASE>.h5' '/hist_lhv_heat' 'wtmean'  0.0  '/lhv_heat_lp_ts'      'b_ZT'   ''      ''        ''       'ge' 0 }
        { 'HDF5/TsAveragedData/hist_lhv_heat_lp_<CASE>.h5' '/hist_lhv_heat' 'wtmean'  0.0  '/lhv_heat_lp'         'b_Zt'   ''      ''        ''       'ge' 0 }
      }
      'DIAGS/hist_meas_ts_lat_heat_<CASE>.h5'
    }


    };

  Nsets = length(MeasSets);

  for icase = 1:Ncases
    Case = CaseList{icase};

    fprintf('*****************************************************************\n');
    fprintf('Generating histogram measurements:\n');
    fprintf('  Case: %s\n', Case);
    fprintf('\n');


    for iset = 1:Nsets
      MeasName     = MeasSets{iset}{1};
      MeasList     = MeasSets{iset}{2};
      OutFtemplate = MeasSets{iset}{3};

      fprintf('    Measurement set: %s\n', MeasName);
      fprintf('\n');

      % Put all measurements into one file per case
      % If the file exists, remove it so that the HDF5 commands
      % can effectively re-create datasets.
      OutFile = regexprep(OutFtemplate, '<CASE>', Case);
      if (exist(OutFile, 'file') == 2)
        delete(OutFile);
      end
  
      Nmeas = length(MeasList);
      for imeas = 1:Nmeas
        Ftemplate = MeasList{imeas}{1};
        Vname     = MeasList{imeas}{2};
        Rmethod   = MeasList{imeas}{3};
        Param     = MeasList{imeas}{4};
        OutVname  = MeasList{imeas}{5};
        Rspec     = num2cell(MeasList{imeas}{6});
        Rrange    = MeasList{imeas}{7};
        Trange    = MeasList{imeas}{8};
        Zrange    = MeasList{imeas}{9};
        SelectOp  = MeasList{imeas}{10};
        SelectVal = MeasList{imeas}{11};

        InFile = regexprep(Ftemplate, '<CASE>', Case);

        % Parse the reduction spec. Always has 4 characters that correspond to x, y, z, t.
        Rindex = 0;  % if index remains zero, then corresponding quantity is not in the input
        Bindex = 0;
        Zindex = 0;
        Tindex = 0;

        Rreduce = false;
        Breduce = false;
        Zreduce = false;
        Treduce = false;

        InForm = {};
        Ndims = 0;
        InDims = { 'x' 'y' 'z' 't' }';
        for i = 1:length(Rspec)
          switch(Rspec{i})
            case { 'r' 'R' }
              RinVar = InDims{i};
              Ndims = Ndims + 1;
              Rindex = Ndims;
              InForm{Ndims} = 'r';
              if (strcmp(Rspec{i}, 'r'))
                Rreduce = true;
              end

            case { 'b' 'B' }
              BinVar = InDims{i};
              Ndims = Ndims + 1;
              Bindex = Ndims;
              InForm{Ndims} = 'b';
              if (strcmp(Rspec{i}, 'b'))
                Breduce = true;
              end

            case { 'z' 'Z' }
              ZinVar = InDims{i};
              Ndims = Ndims + 1;
              Zindex = Ndims;
              InForm{Ndims} = 'z';
              if (strcmp(Rspec{i}, 'z'))
                Zreduce = true;
              end

            case { 't' 'T' }
              TinVar = InDims{i};
              Ndims = Ndims + 1;
              Tindex = Ndims;
              InForm{Ndims} = 't';
              if (strcmp(Rspec{i}, 't'))
                Treduce = true;
              end
          end
        end

        fprintf('      Reading: %s (%s)\n', InFile, Vname);
        fprintf('        Input form: (');
        fprintf('%s', InForm{1});
        for i = 2:Ndims
          fprintf(',%s', InForm{i});
        end
        fprintf(')\n');
        fprintf('\n');

        HDATA = squeeze(h5read(InFile, Vname));
        X     = squeeze(h5read(InFile, '/x_coords'));
        Y     = squeeze(h5read(InFile, '/y_coords'));
        Z     = squeeze(h5read(InFile, '/z_coords'));
        T     = squeeze(h5read(InFile, '/t_coords'));
  
        % Determine indices corresponding to r,b,z,t range specs
        % Default is entire range of dimension sizes
        clear R1;
        clear R2;
        clear B1;
        clear B2;
        clear Z1
        clear Z2
        clear T1;
        clear T2;

        fprintf('        Selection:\n');
        % RADIUS
        if (Rindex > 0)
          % convert to km
          switch(RinVar)
            case 'x'
              R = X ./ 1000;
            case 'y'
              R = Y ./ 1000;
            case 'z'
              R = Z ./ 1000;
            case 't'
              R = T ./ 1000;
          end
          R1 = 1;
          R2 = length(R);
          switch(Rrange)
            case 'core'
              R1 = find(R >= RstartCore, 1, 'first');
              R2 = find(R <= RendCore,   1, 'last');

            case 'rband'
              R1 = find(R >= RstartRband, 1, 'first');
              R2 = find(R <= RendRband,   1, 'last');
  
            case 'env'
              R1 = find(R >= RstartEnv, 1, 'first');
              R2 = find(R <= RendEnv,   1, 'last');
          end

          fprintf('          Radius: "%s" -> %.2f to %.2f (%d:%d)\n', Rrange, R(R1), R(R2), R1, R2);

          % trim the radius coordinates
          R  = R(R1:R2);
          Nr = length(R);

          % add to selection spec
          if (Rindex == 1)
            SelectSpec = sprintf('%d:%d', R1, R2);
          else
            SelectSpec = sprintf('%s,%d:%d', SelectSpec, R1, R2);
          end
        end

        % BINS
        if (Bindex > 0)
          switch(BinVar)
            case 'x'
              B = X;
            case 'y'
              B = Y;
            case 'z'
              B = Z;
            case 't'
              B = T;
          end
          B1 = 1;         
          B2 = length(B);
          switch(SelectOp)
            case 'ge'
              B1 = find(B >= SelectVal, 1, 'first');
  
            case 'le'
              B2 = find(B <= SelectVal, 1, 'last');
          end

          fprintf('          Bin: "%s %.2f" -> %.2f to %.2f (%d:%d)\n', SelectOp, SelectVal, B(B1), B(B2), B1, B2);

          % trim the bins coordinates
          B  = B(B1:B2);
          Nb = length(B);

          % add to selection spec
          if (Bindex == 1)
            SelectSpec = sprintf('%d:%d', B1, B2);
          else
            SelectSpec = sprintf('%s,%d:%d', SelectSpec, B1, B2);
          end
        end

        % HEIGHT
        if (Zindex > 0)
          % convert to km
          switch(ZinVar)
            case 'x'
              H = X ./ 1000;
            case 'y'
              H = Y ./ 1000;
            case 'z'
              H = Z ./ 1000;
            case 't'
              H = T ./ 1000;
          end
          Z1 = 1;
          Z2 = length(H);
          switch(Zrange)
            case { 'maxlev' 'minlev' 'sfc' }
              Z1 = find(H >= Zsfc, 1, 'first');
              Z2 = find(H <= Ztop, 1, 'last');
          end

          fprintf('          Height: "%s" -> %.2f to %.2f (%d:%d)\n', Zrange, H(Z1), H(Z2), Z1, Z2);

          % trim the height coordinates
          H  = H(Z1:Z2);
          Nz = length(H);

          % add to selection spec
          if (Zindex == 1)
            SelectSpec = sprintf('%d:%d', Z1, Z2);
          else
            SelectSpec = sprintf('%s,%d:%d', SelectSpec, Z1, Z2);
          end
        end

        % SIM TIME
        if (Tindex > 0)
          % convert to hours starting with zero
          switch(TinVar)
            case 'x'
              ST = (X ./ 3600) - 42;
            case 'y'
              ST = (Y ./ 3600) - 42;
            case 'z'
              ST = (Z ./ 3600) - 42;
            case 't'
              ST = (T ./ 3600) - 42;
          end
          T1 = 1;
          T2 = length(ST);
          switch(Trange)
            case 'pre_sal'
              T1 = find(ST >= TstartPreSal, 1, 'first');
              T2 = find(ST <= TendPreSal,   1, 'last');
  
            case 'sal'
              T1 = find(ST >= TstartSal, 1, 'first');
              T2 = find(ST <= TendSal,   1, 'last');
          end

          fprintf('          Time: "%s" -> %.2f to %.2f (%d:%d)\n', Trange, ST(T1), ST(T2), T1, T2);

          % trim the sim time coordinates
          ST  = ST(T1:T2);
          Nt = length(ST);

          % add to selection spec
          if (Tindex == 1)
            SelectSpec = sprintf('%d:%d', T1, T2);
          else
            SelectSpec = sprintf('%s,%d:%d', SelectSpec, T1, T2);
          end
        end

        % Trim down the input data according to the selection indices
        SelectCmd = sprintf('HDATA(%s)', SelectSpec);
        MEAS = eval(SelectCmd);
        fprintf('\n');
  
        % If first measurement, then write out coordinates for later
        % use in attaching vars to them. Write out the original coordinates.
        % The original coordinates with their original lengths will always be the
        % appropriate vectors to use for attaching to variables.
        if (imeas == 1)
          Xname = '/x_coords';
          Yname = '/y_coords';
          Zname = '/z_coords';
          Tname = '/t_coords';
  
          CreateDimensionsXyzt(OutFile, X, Y, Z, T, Xname, Yname, Zname, Tname);
          % Add COARDS annotations
          NotateDimensionsXyzt(OutFile, Xname, Yname, Zname, Tname);
        end

        % Reduce the input data. Go in this order:
        %    1. Sum histogram counts (dims r and t)
        %    2. Reduce histogram counts (dim b)
        %    3. Select z level (dim z)
        %
        % Keep original dimensions intact until finished so that [RBZT]index vars
        % can remain constant.
        fprintf('        Reduction:\n');

        if (Rreduce)
          fprintf('          Summing radial counts\n');
          MEAS = squeeze(sum(MEAS, Rindex));

          % Adjust indices, number of dimensions
          if (Bindex > Rindex)
            Bindex = Bindex - 1;
          end
          if (Zindex > Rindex)
            Zindex = Zindex - 1;
          end
          if (Tindex > Rindex)
            Tindex = Tindex - 1;
          end
          Rindex = 0;
          Ndims = Ndims - 1;
        end

        if (Treduce)
          fprintf('          Summing temporal counts\n');
          MEAS = squeeze(sum(MEAS, Tindex));

          % Adjust indices, number of dimensions
          if (Rindex > Tindex)
            Rindex = Rindex - 1;
          end
          if (Bindex > Tindex)
            Bindex = Bindex - 1;
          end
          if (Zindex > Tindex)
            Zindex = Zindex - 1;
          end
          Tindex = 0;
          Ndims = Ndims - 1;
        end

        if (Breduce)
          fprintf('          Reducing bin counts\n');
          MEAS = squeeze(ReduceHists(MEAS, Bindex, B, Rmethod, Param));

          % Adjust indices, number of dimensions
          if (Rindex > Bindex)
            Rindex = Rindex - 1;
          end
          if (Zindex > Bindex)
            Zindex = Zindex - 1;
          end
          if (Tindex > Bindex)
            Tindex = Tindex - 1;
          end
          Bindex = 0;
          Ndims = Ndims - 1;
        end

        if (Zreduce)
          fprintf('          Reducing height\n');
          % Z has been trimmed down to the range we are restricting the search
          % for min or max values, therefore z == 1 represents the surface.
          % Determine which level needs to be selected.
          switch(Zrange)
            case 'sfc'
              ZSEL = 1;

            case 'minlev'
              % 1. grab linear index into MEAS where min occurs
              % 2. translate linear index back to 4D indices, (r,b,z,t)
              %    doesn't matter if MEAS is less than 4D, the extra dims get 1's
              % 3. pick off the index that corresponds to the z dimension
              [ MVAL MLOC ] = min(MEAS(:));
              [ I1 I2 I3 I4 ] = ind2sub(size(MEAS),MLOC);
              MIND = [ I1 I2 I3 I4 ];
              ZSEL = MIND(Zindex);

            case 'maxlev'
              % 1. grab linear index into MEAS where max occurs
              % 2. translate linear index back to 4D indices, (r,b,z,t)
              %    doesn't matter if MEAS is less than 4D, the extra dims get 1's
              % 3. pick off the index that corresponds to the z dimension
              [ MVAL MLOC ] = max(MEAS(:));
              [ I1 I2 I3 I4 ] = ind2sub(size(MEAS),MLOC);
              MIND = [ I1 I2 I3 I4 ];
              ZSEL = MIND(Zindex);
          end
          fprintf('            Z level selected for reduction %s: %.2f (%d)\n', Zrange, H(ZSEL), ZSEL);

          % select the single z level, thus reducing the z dimension
          for i = 1:Ndims
            if (i == Zindex)
              DimSpecs{i} = sprintf('%d', ZSEL);
            else
              DimSpecs{i} = ':';
            end
          end
          SelectCmd = sprintf('squeeze(MEAS(%s))', strjoin(DimSpecs, ','));
          MEAS = eval(SelectCmd);

          % Adjust indices, number of dimensions
          if (Rindex > Zindex)
            Rindex = Rindex - 1;
          end
          if (Bindex > Zindex)
            Bindex = Bindex - 1;
          end
          if (Tindex > Zindex)
            Tindex = Tindex - 1;
          end
          Zindex = 0;
          Ndims = Ndims - 1;
        end
        fprintf('\n');

        % Replace nans with zeros - help make profiles look better
        MEAS(isnan(MEAS)) = 0;

        % Write out measurement
        fprintf('      Writing: %s (%s)\n', OutFile, OutVname)
        fprintf('\n');
  
        % Figure out output dimensions, sizes
        clear OutSize;
        clear DimOrder;

        Oindex = 0;

        for i = 1:length(Rspec)
          switch(Rspec{i})
            case('R')
              Oindex = Oindex + 1;
              OutSize(Oindex) = Nr;
              DimOrder{Oindex} = RinVar;

            case('B')
              Oindex = Oindex + 1;
              OutSize(Oindex) = Nb;
              DimOrder{Oindex} = BinVar;

            case('Z')
              Oindex = Oindex + 1;
              OutSize(Oindex) = Nz;
              DimOrder{Oindex} = ZinVar;

            case('T')
              Oindex = Oindex + 1;
              OutSize(Oindex) = Nt;
              DimOrder{Oindex} = TinVar;

          end
        end

        h5create(OutFile, OutVname, OutSize);
        h5write(OutFile, OutVname, MEAS);

        % Attach dimensions
        AttachDimensionsXyzt(OutFile, OutVname, DimOrder, Xname, Yname, Zname, Tname);
      end % measurements
    end % sets
  end % cases
end % function
