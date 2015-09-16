function [ ] = GenTsdHistMeas()

  Adir = 'AzAveragedData';
  Ddir = 'DIAGS';

  % make sure output directory exists
  if (exist(Ddir, 'dir') ~= 7)
    mkdir(Ddir);
  end

  % list of simulation cases
  CaseList = {
    'TSD_SAL_DUST'
    'TSD_SAL_NODUST'
    'TSD_NONSAL_DUST'
    'TSD_NONSAL_NODUST'
    };
  Ncases = length(CaseList);

  % Description of measurements
  %   <measure_name> <measure_list> <out_file_prefix>
  %
  %   where <measure_list> is one or more of:
  %     <file_prefix> <var_name> <reduction_method> <arg_for_reduction_method> <out_var_name> <select_op> <select_val>
  MeasSets = {
%    % storm speed measurements
%    {
%      'Speed'
%      {
%        { 'hist_speed' '/speed' 'farea'  0.50 '/avg_speed_fa'     'ge' 0   }
%        { 'hist_speed' '/speed' 'wtmean'  0.0 '/avg_speed_wm'     'ge' 0   }
%        { 'hist_speed' '/speed' 'farea'  0.99 '/max_speed_fa'     'ge' 0   }
%  
%        { 'hist_speed10m' '/speed10m' 'farea'  0.50 '/avg_speed10m_fa'     'ge' 0   }
%        { 'hist_speed10m' '/speed10m' 'wtmean'  0.0 '/avg_speed10m_wm'     'ge' 0   }
%        { 'hist_speed10m' '/speed10m' 'farea'  0.99 '/max_speed10m_fa'     'ge' 0   }
%  
%        { 'hist_speed_t' '/speed_t' 'farea'  0.50 '/avg_speed_t_fa'     'ge' -100   }
%        { 'hist_speed_t' '/speed_t' 'wtmean'  0.0 '/avg_speed_t_wm'     'ge' -100   }
%        { 'hist_speed_t' '/speed_t' 'farea'  0.99 '/max_speed_t_fa'     'ge' -100   }
%  
%        { 'hist_speed_r' '/speed_r' 'farea'  0.50 '/avg_speed_r_fa'     'ge' -100   }
%        { 'hist_speed_r' '/speed_r' 'wtmean'  0.0 '/avg_speed_r_wm'     'ge' -100   }
%        { 'hist_speed_r' '/speed_r' 'farea'  0.99 '/max_speed_r_fa'     'ge' -100   }
%      }
%      'hist_meas_speed'
%    }
%  
%    % storm pressure measurements
%    {
%      'Pressure'
%      {
%        { 'hist_press' '/press' 'farea'  0.50 '/avg_press_fa'     'ge' 0   }
%        { 'hist_press' '/press' 'wtmean'  0.0 '/avg_press_wm'     'ge' 0   }
%        { 'hist_press' '/press' 'farea'  0.01 '/min_press_fa'     'ge' 0   }
%  
%        { 'hist_sea_press' '/sea_press' 'farea'  0.50 '/avg_sea_press_fa'     'ge' 0   }
%        { 'hist_sea_press' '/sea_press' 'wtmean'  0.0 '/avg_sea_press_wm'     'ge' 0   }
%        { 'hist_sea_press' '/sea_press' 'farea'  0.01 '/min_sea_press_fa'     'ge' 0   }
%      }
%      'hist_meas_press'
%    }
%
%    % precip rate measurements
%    {
%      'Precip Rate'
%      {
%        { 'hist_pcprate' '/pcprate' 'farea'  0.50 '/avg_pcprate_fa'     'ge' 0   }
%        { 'hist_pcprate' '/pcprate' 'wtmean'  0.0 '/avg_pcprate_wm'     'ge' 0   }
%        { 'hist_pcprate' '/pcprate' 'farea'  0.99 '/max_pcprate_fa'     'ge' 0   }
%      }
%      'hist_meas_pcprate'
%    }
%
%    % vertially integrated condensate measurements
%    {
%      'Vert Cond'
%      {
%        { 'hist_vint_cond' '/vint_cond' 'farea'  0.50 '/avg_vint_cond_fa'     'ge' 0   }
%        { 'hist_vint_cond' '/vint_cond' 'wtmean'  0.0 '/avg_vint_cond_wm'     'ge' 0   }
%        { 'hist_vint_cond' '/vint_cond' 'farea'  0.99 '/max_vint_cond_fa'     'ge' 0   }
%      }
%      'hist_meas_vint_cond'
%    }
%
%    % vertical velocity measurements
%    {
%      'Vertical Velocity'
%      {
%        { 'hist_up' '/w' 'farea'  0.50 '/avg_updraft_fa'     'ge'  0.01   }
%        { 'hist_up' '/w' 'wtmean'  0.0 '/avg_updraft_wm'     'ge'  0.01   }
%        { 'hist_w'  '/w' 'farea'  0.99 '/max_updraft_fa'     'ge'  0.05   }
%
%        { 'hist_dn' '/w' 'farea'  0.50 '/avg_dndraft_fa'     'le' -0.01   }
%        { 'hist_dn' '/w' 'wtmean'  0.0 '/avg_dndraft_wm'     'le' -0.01   }
%        { 'hist_w'  '/w' 'farea'  0.01 '/max_dndraft_fa'     'le' -0.05   }
%      }
%      'hist_meas_w'
%    }
%
%    % theta_e measurements
%    {
%      'Theta-E'
%      {
%        { 'hist_theta_e' '/theta_e' 'farea'  0.50 '/avg_theta_e_fa'     'ge' 0   }
%        { 'hist_theta_e' '/theta_e' 'wtmean'  0.0 '/avg_theta_e_wm'     'ge' 0   }
%        { 'hist_theta_e' '/theta_e' 'farea'  0.99 '/max_theta_e_fa'     'ge' 0   }
%      }
%      'hist_meas_theta_e'
%    }
%
%    % dust measurements
%    {
%      'Dust'
%      {
%        { 'hist_d1_mass' '/d1_mass' 'farea'  0.50 '/avg_d1_mass_fa'     'ge' 0   }
%        { 'hist_d1_mass' '/d1_mass' 'wtmean'  0.0 '/avg_d1_mass_wm'     'ge' 0   }
%        { 'hist_d1_mass' '/d1_mass' 'farea'  0.99 '/max_d1_mass_fa'     'ge' 0   }
%
%        { 'hist_d1_num' '/d1_num' 'farea'  0.50 '/avg_d1_num_fa'     'ge' 0   }
%        { 'hist_d1_num' '/d1_num' 'wtmean'  0.0 '/avg_d1_num_wm'     'ge' 0   }
%        { 'hist_d1_num' '/d1_num' 'farea'  0.99 '/max_d1_num_fa'     'ge' 0   }
%
%        { 'hist_d2_mass' '/d2_mass' 'farea'  0.50 '/avg_d2_mass_fa'     'ge' 0   }
%        { 'hist_d2_mass' '/d2_mass' 'wtmean'  0.0 '/avg_d2_mass_wm'     'ge' 0   }
%        { 'hist_d2_mass' '/d2_mass' 'farea'  0.99 '/max_d2_mass_fa'     'ge' 0   }
%
%        { 'hist_d2_num' '/d2_num' 'farea'  0.50 '/avg_d2_num_fa'     'ge' 0   }
%        { 'hist_d2_num' '/d2_num' 'wtmean'  0.0 '/avg_d2_num_wm'     'ge' 0   }
%        { 'hist_d2_num' '/d2_num' 'farea'  0.99 '/max_d2_num_fa'     'ge' 0   }
%      }
%      'hist_meas_dust'
%    }

    % cloud
    {
      'Cloud'
     {
        { 'hist_cloud' '/cloud' 'farea'  0.50 '/avg_cloud_fa'     'ge' 0   }
        { 'hist_cloud' '/cloud' 'wtmean'  0.0 '/avg_cloud_wm'     'ge' 0   }
        { 'hist_cloud' '/cloud' 'farea'  0.99 '/max_cloud_fa'     'ge' 0   }
      }
      'hist_meas_cloud'
    }

    % rain
    {
      'Rain'
     {
        { 'hist_rain' '/rain' 'farea'  0.50 '/avg_rain_fa'     'ge' 0   }
        { 'hist_rain' '/rain' 'wtmean'  0.0 '/avg_rain_wm'     'ge' 0   }
        { 'hist_rain' '/rain' 'farea'  0.99 '/max_rain_fa'     'ge' 0   }
      }
      'hist_meas_rain'
    }

    % pristine
    {
      'Pristine'
     {
        { 'hist_pris' '/pris' 'farea'  0.50 '/avg_pris_fa'     'ge' 0   }
        { 'hist_pris' '/pris' 'wtmean'  0.0 '/avg_pris_wm'     'ge' 0   }
        { 'hist_pris' '/pris' 'farea'  0.99 '/max_pris_fa'     'ge' 0   }
      }
      'hist_meas_pris'
    }

    % snow
    {
      'Snow'
     {
        { 'hist_snow' '/snow' 'farea'  0.50 '/avg_snow_fa'     'ge' 0   }
        { 'hist_snow' '/snow' 'wtmean'  0.0 '/avg_snow_wm'     'ge' 0   }
        { 'hist_snow' '/snow' 'farea'  0.99 '/max_snow_fa'     'ge' 0   }
      }
      'hist_meas_snow'
    }

    % aggregates
    {
      'Aggregates'
     {
        { 'hist_aggr' '/aggr' 'farea'  0.50 '/avg_aggr_fa'     'ge' 0   }
        { 'hist_aggr' '/aggr' 'wtmean'  0.0 '/avg_aggr_wm'     'ge' 0   }
        { 'hist_aggr' '/aggr' 'farea'  0.99 '/max_aggr_fa'     'ge' 0   }
      }
      'hist_meas_aggr'
    }

    % graupel
    {
      'Graupel'
     {
        { 'hist_graup' '/graup' 'farea'  0.50 '/avg_graup_fa'     'ge' 0   }
        { 'hist_graup' '/graup' 'wtmean'  0.0 '/avg_graup_wm'     'ge' 0   }
        { 'hist_graup' '/graup' 'farea'  0.99 '/max_graup_fa'     'ge' 0   }
      }
      'hist_meas_graup'
    }

    % hail
    {
      'Hail'
     {
        { 'hist_hail' '/hail' 'farea'  0.50 '/avg_hail_fa'     'ge' 0   }
        { 'hist_hail' '/hail' 'wtmean'  0.0 '/avg_hail_wm'     'ge' 0   }
        { 'hist_hail' '/hail' 'farea'  0.99 '/max_hail_fa'     'ge' 0   }
      }
      'hist_meas_hail'
    }

%    % cooling via latent heat of freezing
%    {
%      'LHF Cooling'
%     {
%        { 'hist_lhf_cool' '/lhf_cool' 'farea'  0.50 '/avg_lhf_cool_fa'     'le' 0   }
%        { 'hist_lhf_cool' '/lhf_cool' 'wtmean'  0.0 '/avg_lhf_cool_wm'     'le' 0   }
%        { 'hist_lhf_cool' '/lhf_cool' 'farea'  0.99 '/max_lhf_cool_fa'     'le' 0   }
%      }
%      'hist_meas_lhf_cool'
%    }
%
%    % heating via latent heat of freezing
%    {
%      'LHF Heating'
%     {
%        { 'hist_lhf_heat' '/lhf_heat' 'farea'  0.50 '/avg_lhf_heat_fa'     'ge' 0   }
%        { 'hist_lhf_heat' '/lhf_heat' 'wtmean'  0.0 '/avg_lhf_heat_wm'     'ge' 0   }
%        { 'hist_lhf_heat' '/lhf_heat' 'farea'  0.99 '/max_lhf_heat_fa'     'ge' 0   }
%      }
%      'hist_meas_lhf_heat'
%    }
%
%    % cooling via latent heat of vaporization
%    {
%      'LHV Cooling'
%     {
%        { 'hist_lhv_cool' '/lhv_cool' 'farea'  0.50 '/avg_lhv_cool_fa'     'le' 0   }
%        { 'hist_lhv_cool' '/lhv_cool' 'wtmean'  0.0 '/avg_lhv_cool_wm'     'le' 0   }
%        { 'hist_lhv_cool' '/lhv_cool' 'farea'  0.99 '/max_lhv_cool_fa'     'le' 0   }
%      }
%      'hist_meas_lhv_cool'
%    }
%
%    % heating via latent heat of vaporization
%    {
%      'LHV Heating'
%     {
%        { 'hist_lhv_heat' '/lhv_heat' 'farea'  0.50 '/avg_lhv_heat_fa'     'ge' 0   }
%        { 'hist_lhv_heat' '/lhv_heat' 'wtmean'  0.0 '/avg_lhv_heat_wm'     'ge' 0   }
%        { 'hist_lhv_heat' '/lhv_heat' 'farea'  0.99 '/max_lhv_heat_fa'     'ge' 0   }
%      }
%      'hist_meas_lhv_heat'
%    }

    };

  Nsets = length(MeasSets);

  for icase = 1:Ncases
    Case = CaseList{icase};

    fprintf('*****************************************************************\n');
    fprintf('Generating histogram measurements:\n');
    fprintf('  Case: %s\n', Case);
    fprintf('\n');


    for iset = 1:Nsets
      MeasName   = MeasSets{iset}{1};
      MeasList   = MeasSets{iset}{2};
      OutFprefix = MeasSets{iset}{3};

      fprintf('    Measurement set: %s\n', MeasName);
      fprintf('\n');

      % Put all measurements into one file per case
      % If the file exists, remove it so that the HDF5 commands
      % can effectively re-create datasets.
      OutFile = sprintf('%s/%s_%s.h5', Ddir, OutFprefix, Case);
      if (exist(OutFile, 'file') == 2)
        delete(OutFile);
      end
  
      Nmeas = length(MeasList);
      MaxNz = 0;
      for imeas = 1:Nmeas
        Fprefix   = MeasList{imeas}{1};
        Vname     = MeasList{imeas}{2};
        Rmethod   = MeasList{imeas}{3};
        Param     = MeasList{imeas}{4};
        OutVname  = MeasList{imeas}{5};
        SelectOp  = MeasList{imeas}{6};
        SelectVal = MeasList{imeas}{7};
  
        InFile = sprintf('%s/%s_%s.h5', Adir, Fprefix, Case);
  
        fprintf('      Reading: %s (%s)\n', InFile, Vname);
        fprintf('        Reduction method: %s (%.2f)\n', Rmethod, Param);
        fprintf('        Selection: %s %.2f\n', SelectOp, SelectVal);
  
        % Read in data which will be 4D -> (x,y,z,t)
        %
        %     x --> radial bands
        %     y --> histogram bins
        %     z --> height
        %     t --> time
        %
        HDATA = squeeze(h5read(InFile, Vname));
        BINS  = squeeze(h5read(InFile, '/y_coords'));
  
        % Assume same r,z,t values for all measurements
        X    = squeeze(h5read(InFile, '/x_coords'));
        Y    = 1; % dummy dimension for output
        Z    = squeeze(h5read(InFile, '/z_coords'));
        T    = squeeze(h5read(InFile, '/t_coords'));
  
        Nx = length(X);
        Ny = length(Y);
        Nz = length(Z);
        Nt = length(T);

        % Keep track of the coordinates of the full (maximum sized) z dimension
        if (Nz > MaxNz)
          MaxNz = Nz;
          AllNz = Nz;
          AllZ = Z;
        end
  
        % Reduce the histograms level by level to preserve vertical structure
        % MEAS will be organized as: (r,z,t)
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
  
        MEAS = squeeze(ReduceHists(HDATA(:,B1:B2,:,:), 2, BINS(B1:B2), Rmethod, Param));
  
        % Write out measurement
        fprintf('      Writing: %s (%s)\n', OutFile, OutVname)
        fprintf('\n');
  
        % Write out measurement -> force to be (x,y,z,t) for dimension
        % attach code below.
        OutVar = reshape(MEAS, [ Nx Ny Nz Nt ]);
        h5create(OutFile, OutVname, size(OutVar));
        h5write(OutFile, OutVname, OutVar);
      end % measurements
  
      % Create the dimensions
      fprintf('      Creating dimensions: %s\n', OutFile);
      fprintf('\n');
  
      Xname = '/x_coords';
      Yname = '/y_coords';
      Zname = '/z_coords';
      Tname = '/t_coords';
  
      CreateDimensionsXyzt(OutFile, X, Y, AllZ, T, Xname, Yname, Zname, Tname);
  
      % Attach dimensions to all variables
      DimOrder = { 'x' 'y' 'z' 't' };
      for imeas = 1:Nmeas
        Vname = MeasList{imeas}{5};   % use the output var name
        AttachDimensionsXyzt(OutFile, Vname, DimOrder, Xname, Yname, Zname, Tname);
      end

      % Add COARDS annotations
      NotateDimensionsXyzt(OutFile, Xname, Yname, Zname, Tname);

    end % sets
  end % cases
end % function
