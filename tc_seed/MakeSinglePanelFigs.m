function [] = MakeSinglePanelFigs(ConfigFile)
% MakeSinglePanelFigs - covert matlab fig to JPEG

% Read in the config data
[ Config ] = ReadConfig(ConfigFile);

PlotDir = Config.PlotDir;
FigDir = Config.FigDir;

% make sure output directory exists
if (exist(FigDir, 'dir') ~= 7)
  mkdir(FigDir);
end

% list of drawings to translate

%FigList = { 'InitVortex'
%            'SampleCcnProf'
%            'SpinUpVtSpl'
%            'StormPhases'
%            'StormRegions'
%            'KeVt'
%            'TS_rmw' };

FigList = { 'EOF1_theta_e_TCS_SD_C0500'
            'EOF1_theta_e_TCS_SD_C2000'
            'PC1_theta_e_TCS_SD_C0500'
            'PC1_theta_e_TCS_SD_C2000' };

%FigList = { 'prs_vaptott_vt0p5_EW_RI_TCS_GN_C0500'
%            'prs_vaptott_vt0p5_EW_RI_TCS_GN_C1000'
%            'prs_vaptott_vt0p5_EW_RI_TCS_GN_C2000'
%            'prs_lh_tott_lht1p0_EW_RI_TCS_GN_C0500'
%            'prs_lh_tott_lht1p0_EW_RI_TCS_GN_C1000'
%            'prs_lh_tott_lht1p0_EW_RI_TCS_GN_C2000'
%            'prs_w_EW_RI_TCS_GN_C0500'
%            'prs_w_EW_RI_TCS_GN_C1000'
%            'prs_w_EW_RI_TCS_GN_C2000'
%            'prs_vaptott_vt0p5_EW_SS_TCS_GN_C0500'
%            'prs_vaptott_vt0p5_EW_SS_TCS_GN_C1000'
%            'prs_vaptott_vt0p5_EW_SS_TCS_GN_C2000'
%            'prs_lh_tott_lht1p0_EW_SS_TCS_GN_C0500'
%            'prs_lh_tott_lht1p0_EW_SS_TCS_GN_C1000'
%            'prs_lh_tott_lht1p0_EW_SS_TCS_GN_C2000'
%            'prs_w_EW_SS_TCS_GN_C0500'
%            'prs_w_EW_SS_TCS_GN_C1000'
%            'prs_w_EW_SS_TCS_GN_C2000'
%            'prs_relvortz_rvz1em4_EW_RI_TCS_GN_C0100'
%            'prs_relvortz_rvz1em4_EW_RI_TCS_GN_C0500'
%            'prs_relvortz_rvz1em4_EW_RI_TCS_GN_C1000'
%            'prs_relvortz_rvz1em4_EW_RI_TCS_GN_C2000'
%            'prs_relvortz_rvz1em4_EW_SS_TCS_GN_C0100'
%            'prs_relvortz_rvz1em4_EW_SS_TCS_GN_C0500'
%            'prs_relvortz_rvz1em4_EW_SS_TCS_GN_C1000'
%            'prs_relvortz_rvz1em4_EW_SS_TCS_GN_C2000'
%            'prs_relvortz_rvz1em4_diff_EW_RI_TCS_GN_C0500'
%            'prs_relvortz_rvz1em4_diff_EW_RI_TCS_GN_C1000'
%            'prs_relvortz_rvz1em4_diff_EW_RI_TCS_GN_C2000'
%            'prs_relvortz_rvz1em4_diff_EW_SS_TCS_GN_C0500'
%            'prs_relvortz_rvz1em4_diff_EW_SS_TCS_GN_C1000'
%            'prs_relvortz_rvz1em4_diff_EW_SS_TCS_GN_C2000' };

for i = 1:length(FigList)
  InFile  = sprintf('%s/%s.fig', PlotDir, FigList{i});
  OutFile = sprintf('%s/%s.jpg', FigDir,  FigList{i});

  fprintf('MATLAB figure file: %s\n', InFile);
  fprintf('Writing file: %s\n', OutFile);
  fprintf('\n');

  Fig = openfig(InFile);
  saveas(Fig, OutFile);
  close(Fig);
end

end
