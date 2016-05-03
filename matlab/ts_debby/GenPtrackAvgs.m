function [ ] = GenPtrackAvgs()

  CaseList = {
   'TSD_SAL_DUST'
   'TSD_SAL_NODUST'
   'TSD_NONSAL_DUST'
   'TSD_NONSAL_NODUST'
   };
  Nc = length(CaseList);

  MeasList = {
    { 'XsectionData/ptrack_theta_<CASE>.h5' '/theta' 'pre_sal' '/ps_theta' }
    { 'XsectionData/ptrack_theta_<CASE>.h5' '/theta' 'sal'     '/s_theta'  }

    { 'XsectionData/ptrack_theta_e_<CASE>.h5' '/theta_e' 'pre_sal' '/ps_theta_e' }
    { 'XsectionData/ptrack_theta_e_<CASE>.h5' '/theta_e' 'sal'     '/s_theta_e'  }

    { 'XsectionData/ptrack_tempc_<CASE>.h5' '/tempc' 'pre_sal' '/ps_tempc' }
    { 'XsectionData/ptrack_tempc_<CASE>.h5' '/tempc' 'sal'     '/s_tempc'  }
    };
  Nm = length(MeasList);

  PreSalStart = 10;  % start, end Pre-SAL time period, sim time in hours
  PreSalEnd   = 30;
  SalStart    = 40;  % start, end SAL time period
  SalEnd      = 60;

  Xname = '/x_coords';
  Yname = '/y_coords';
  Zname = '/z_coords';
  Tname = '/t_coords';

  fprintf('Generating ptrack averages\n');
  fprintf('\n');

  for icase = 1:Nc
    Case = CaseList{icase};
    fprintf('  Case: %s\n', Case);
    fprintf('\n');

    % Remove output file if it exists so that new dataset can be written.
    OutFname = sprintf('DIAGS/ptrack_avgs_%s.h5', Case);
    if (exist(OutFname, 'file') == 2)
      delete(OutFname);
    end

    for imeas = 1:Nm
      InFname   = regexprep(MeasList{imeas}{1}, '<CASE>', Case);
      InVname   = MeasList{imeas}{2};
      AvgPeriod = MeasList{imeas}{3};
      OutVname  = MeasList{imeas}{4};

      % if on first measurement, create the coordinate datasets
      if (imeas == 1)
        X = squeeze(h5read(InFname, Xname));
        Y = squeeze(h5read(InFname, Yname));
        Z = squeeze(h5read(InFname, Zname));
        T = squeeze(h5read(InFname, Tname));

        Nx = length(X);
        Ny = length(Y);
        Nz = length(Z);
        Nt = length(T);

        % Write coordinates into output file
        CreateDimensionsXyzt(OutFname, X, Y, Z, T, Xname, Yname, Zname, Tname);
        NotateDimensionsXyzt(OutFname, Xname, Yname, Zname, Tname);

        % Create sim time, to be used for selecting time period for averaging
        %   sim time is in hours, starting with zero
        SIM_T = T ./ 3600 - 42;
      end

      fprintf('    Reading: %s (%s)\n', InFname, InVname);
      PDATA = squeeze(h5read(InFname, InVname));
  
      % find the start and end of the time period for averaging
      switch AvgPeriod
        case 'pre_sal'
           T1 = find(SIM_T >= PreSalStart, 1, 'first');
           T2 = find(SIM_T <= PreSalEnd,   1, 'last');
        case 'sal'
           T1 = find(SIM_T >= SalStart, 1, 'first');
           T2 = find(SIM_T <= SalEnd,   1, 'last');
        otherwise
          fprintf('      WARNING: unknown average period: %s, using entire time range\n', AvgPeriod);
          T1 = 1;
          T2 = Nt;
      end
      fprintf('    Average Period: %s -> %.1f (%d) to %.1f (%d)\n', AvgPeriod, SIM_T(T1), T1, SIM_T(T2), T2);

      PAVG = squeeze(mean(PDATA(:,:,T1:T2),3));

      % Write out average
      Vsize = [ Nx Nz ];
      DimOrder = { 'x' 'z' };

      fprintf('    Writing: %s (%s)\n', OutFname, OutVname);
      h5create(OutFname, OutVname, Vsize)
      h5write(OutFname, OutVname, PAVG);
      AttachDimensionsXyzt(OutFname, OutVname, DimOrder, Xname, Yname, Zname, Tname);
      fprintf('\n');
    end
  end
end
