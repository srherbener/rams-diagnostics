% Script to generate sets of profile plots
%
%   domain/time average of latent heat of vaporization
%     CCN == 50/cc and SST == 298K is the control run
%     generate differences from the control run in each set

clear;

SST = [ 293 298 303 ];
%CCN = [ 50 100 200 400 800 1200 1600 ]; % too cluttered
CCN = [ 100 400 800 1600 ];
%Lcolors = { 'k', 'm', 'b', 'c', 'g', 'y', 'r' };
Lcolors = { 'b', 'g', 'y', 'r' };

CNTL_SST = 298;
CNTL_CCN = 50;

Ns = length(SST);
Nc = length(CCN);

Xlabel = 'Latent Heat Difference (K)';
Zlabel = 'Height (m)';

X = (-0.6:.1:0.6);

% Create the output directory if it doesn't exist
OutDir = 'PLOTS';
if (exist(OutDir,'dir') == 0)
  mkdir(OutDir);
end

% read in the control profile
Hfile = sprintf('DIAGS/lh_vapt_tdavg_ATEX_C%04d_S%03d.h5',CNTL_CCN,CNTL_SST);
fprintf('Reading control HDF5 file: %s\n\n', Hfile);
CNTL_LHV = hdf5read(Hfile, '/lh_vapt')';
Zall = hdf5read(Hfile, '/z_coords')'; % Z should be same in all files

Z1 = 2;
Z2 = length(Zall) - 15; % each level is 100m, this chops off the top 1500m
Z = Zall(Z1:Z2);

% Each set: all CCN levels for a single given SST
for i = 1:Ns
  fprintf('Generating plot for SST = %d\n', SST(i));
  Ptitle = sprintf('Latent Heat of Vaporization: SST = %d', SST(i));

  % collect the latent heat data
  for j = 1:Nc
    Hfile = sprintf('DIAGS/lh_vapt_tdavg_ATEX_C%04d_S%03d.h5',CCN(j),SST(i));
    fprintf('  Reading HDF5 file: %s\n', Hfile);
    LHV_DOMAVG = hdf5read(Hfile, '/lh_vapt');
    LHV(j,:) = squeeze(mean(LHV_DOMAVG,2)); % time average

    % subtract off the control profile
    LHV_DIFF(j,:) = LHV(j,Z1:Z2) - CNTL_LHV(Z1:Z2);

    % generate the legend text
    LegText(j) = { sprintf('CCN: %d/cc', CCN(j)) };
  end
  fprintf('\n');

  % plot it
  OutFile = sprintf('%s/lh_vapt_S%03d.jpg', OutDir, SST(i));
  fprintf('  Saving plot in file: %s\n\n', OutFile);
  PlotProfSet(X, Z, LHV_DIFF, Xlabel, Zlabel, Ptitle, Lcolors, LegText, 'NorthEast', OutFile);
end
