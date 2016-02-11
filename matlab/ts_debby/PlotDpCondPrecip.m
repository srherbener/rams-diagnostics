function [ ] = PlotDpCondPrecip()
% PlotDpCondPrecip 

  Pdir = 'Plots';
  if (exist(Pdir, 'dir') ~= 7)
      mkdir(Pdir);
  end

  Fsize = 13;
  
  % Time series of lower level Mdhy
  InFile = 'DIAGS/total_mass_TSD_SAL_DUST.h5';
  InVname = '/sal_tcond_total_mass_llev';
  fprintf('Reading: %s (%s)\n', InFile, InVname);

  TS_SAL_TCOND = squeeze(h5read(InFile, InVname)); % g
  TS_SAL_TCOND = TS_SAL_TCOND .* 1e-15;  % convert to 10^3 Tg

  T = squeeze(h5read(InFile, '/t_coords'))./3600 - 42; % h

  % Time series of surface dust accumulation
  InFile = 'DIAGS/total_mass_TSD_SAL_DUST.h5';
  InVname = '/sal_accpcp_total_mass';
  fprintf('Reading: %s (%s)\n', InFile, InVname);

  TS_SAL_ACCPCP = squeeze(h5read(InFile, InVname)); % g
  TS_SAL_ACCPCP = TS_SAL_ACCPCP .* 1e-15;  % convert to 10^3 Tg

  % plot, 2 panels
  OutFile = sprintf('%s/DpFigCondPrecip.jpg', Pdir);
  
  Fig = figure;

  % create the Hovmoller of sal_dust_hydro
  Paxes = subplot(2,1,1);
  PlotDpFigTseries(Paxes, T, TS_SAL_TCOND, 'a', 'SAL', 'M_c_o_n_d (10^3 Tg)', Fsize, 0, 1);
  
  Paxes = subplot(2,1,2);
  PlotDpFigTseries(Paxes, T, TS_SAL_ACCPCP, 'b', 'SAL', 'M_p_r_e_c_i_p (10^3 Tg)', Fsize, 1, 1);

  fprintf('Writing: %s\n', OutFile);
  saveas(Fig, OutFile);
  close(Fig);
end
