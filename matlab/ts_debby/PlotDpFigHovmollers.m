function [ ] = PlotDpFigHovmollers()

  Pdir = 'Plots';
  if (exist(Pdir, 'dir') ~= 7)
      mkdir(Pdir);
  end
  
  Fsize = 13;

  CaseList = {
    'TSD_SAL_DUST'
    'TSD_SD'
    'TSD_SD_1G'
    };
  Ncases = length(CaseList);

  for icase = 1:Ncases
    Case = CaseList{icase};

    % Hovmoller data for Md
    InFile = sprintf('DIAGS/storm_hovmollers_%s.h5', Case);
    InVname = '/sal_ar_dust_mass';
    fprintf('Reading: %s (%s)\n', InFile, InVname);
    HOV_DUST = squeeze(h5read(InFile, InVname));
    Z    = squeeze(h5read(InFile, '/z_coords'))./1000;         % convert to km
    if (strcmp(Case, 'TSD_SD_1G'))
      T = squeeze(h5read(InFile, '/t_coords'))./3600 - 12;    % convert to sim time in hours
    else
      T = squeeze(h5read(InFile, '/t_coords'))./3600 - 42;    % convert to sim time in hours
    end
  
    % Hovmoller data for Mdhy
    InFile = sprintf('DIAGS/storm_hovmollers_%s.h5', Case);
    InVname = '/sal_ar_dust_hydro';
    fprintf('Reading: %s (%s)\n', InFile, InVname);
    HOV_DHY = squeeze(h5read(InFile, InVname));
  
    % Hovmoller data for Mdrgn
    InFile = sprintf('DIAGS/storm_hovmollers_%s.h5', Case);
    InVname = '/sal_ar_ra_mass';
    fprintf('Reading: %s (%s)\n', InFile, InVname);
    HOV_DRGN = squeeze(h5read(InFile, InVname));
  
    % plot
    Fig = figure;
  
    % Md
    Paxes = subplot(3,1,1);
    PlotDpFigDustHov(Paxes, T, Z, HOV_DUST, 'a', 'M_D', Fsize, 0, 1, 0, 2);
  
    % Mdhy
    Paxes = subplot(3,1,2);
    PlotDpFigDustHov(Paxes, T, Z, HOV_DHY, 'b', 'M_D_H_Y', Fsize, 0, 1, 0, 0);
  
    % Mdrgn
    Paxes = subplot(3,1,3);
    PlotDpFigDustHov(Paxes, T, Z, HOV_DRGN, 'c', 'M_D_R_G_N', Fsize, 1, 1, 0, 1);
    
    OutFile = sprintf('%s/DpFigHovmollers_%s.jpg', Pdir, Case);
    fprintf('Writing: %s\n', OutFile);
    saveas(Fig, OutFile);
    close(Fig);

    fprintf('\n');
  end
end
