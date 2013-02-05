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

%FigList = {
%  'InitVortex'
%  'SampleCcnProf'
%  'SpinUpVtSpl'
%  'StormPhases'
%  'StormRegions'
%  'KeVt'
%  'KeVt_ext'
%  'TS_rmw',
%  'TS_Vint_Vapt_SC'
%  'TS_Vint_Vapt_RB'
%  'TS_Vint_Lht_SC'
%  'TS_Vint_Lht_RB'
%  'TS_Vint_VaptAll_RB'
%  'TS_Vint_LhtAll_RB'
%  'TS_DeltapAll_SC_RB'
%  'vint_cond_TCS_SD_C0100'
%  'vint_cond_TCS_SD_C0500'
%  'vint_cond_TCS_SD_C1000'
%  'vint_cond_TCS_SD_C2000'
%  };

%FigList = {
%  'EOF1_theta_e_TCS_SD_C0500'
%  'EOF1_theta_e_TCS_SD_C1000'
%  'EOF1_theta_e_TCS_SD_C2000'
%  'PC1_theta_e_TCS_SD_C0500'
%  'PC1_theta_e_TCS_SD_C1000'
%  'PC1_theta_e_TCS_SD_C2000'
%  };

%FigList = {
%  'pts_cloud_cm3_twp4_SC_LC_TCS_SD_C0100'
%  'pts_cloud_cm3_twp4_SC_LC_TCS_SD_C0500'
%  'pts_cloud_cm3_twp4_SC_LC_TCS_SD_C1000'
%  'pts_cloud_cm3_twp4_SC_LC_TCS_SD_C2000'
%  'pts_cloud_d_twp4_SC_LC_TCS_SD_C0100'
%  'pts_cloud_d_twp4_SC_LC_TCS_SD_C0500'
%  'pts_cloud_d_twp4_SC_LC_TCS_SD_C1000'
%  'pts_cloud_d_twp4_SC_LC_TCS_SD_C2000'
%  'pts_lh_tott_twp4_lht1p0_SC_LC_TCS_SD_C0100'
%  'pts_lh_tott_twp4_lht1p0_SC_LC_TCS_SD_C0500'
%  'pts_lh_tott_twp4_lht1p0_SC_LC_TCS_SD_C1000'
%  'pts_lh_tott_twp4_lht1p0_SC_LC_TCS_SD_C2000'
%  'pts_cloud_cm3_twp4_RB_LC_TCS_SD_C0100'
%  'pts_cloud_cm3_twp4_RB_LC_TCS_SD_C0500'
%  'pts_cloud_cm3_twp4_RB_LC_TCS_SD_C1000'
%  'pts_cloud_cm3_twp4_RB_LC_TCS_SD_C2000'
%  'pts_cloud_d_twp4_RB_LC_TCS_SD_C0100'
%  'pts_cloud_d_twp4_RB_LC_TCS_SD_C0500'
%  'pts_cloud_d_twp4_RB_LC_TCS_SD_C1000'
%  'pts_cloud_d_twp4_RB_LC_TCS_SD_C2000'
%  'pts_lh_tott_twp4_lht1p0_RB_LC_TCS_SD_C0100'
%  'pts_lh_tott_twp4_lht1p0_RB_LC_TCS_SD_C0500'
%  'pts_lh_tott_twp4_lht1p0_RB_LC_TCS_SD_C1000'
%  'pts_lh_tott_twp4_lht1p0_RB_LC_TCS_SD_C2000'
%  'pts_relhum_RB_LC_TCS_SD_C0100'
%  'pts_relhum_RB_LC_TCS_SD_C0500'
%  'pts_relhum_RB_LC_TCS_SD_C1000'
%  'pts_relhum_RB_LC_TCS_SD_C2000'
%  'pts_theta_e_RB_LC_TCS_SD_C0100'
%  'pts_theta_e_RB_LC_TCS_SD_C0500'
%  'pts_theta_e_RB_LC_TCS_SD_C1000'
%  'pts_theta_e_RB_LC_TCS_SD_C2000'
%  };

%FigList = {
%  'prs_vaptott_vt0p5_SC_RI_TCS_SD_C0500'
%  'prs_vaptott_vt0p5_SC_RI_TCS_SD_C1000'
%  'prs_vaptott_vt0p5_SC_RI_TCS_SD_C2000'
%  'prs_vaptott_vt0p5_SC_TR_TCS_SD_C0500'
%  'prs_vaptott_vt0p5_SC_TR_TCS_SD_C1000'
%  'prs_vaptott_vt0p5_SC_TR_TCS_SD_C2000'
%  'prs_vaptott_vt0p5_SC_SS_TCS_SD_C0500'
%  'prs_vaptott_vt0p5_SC_SS_TCS_SD_C1000'
%  'prs_vaptott_vt0p5_SC_SS_TCS_SD_C2000'
%  'prs_lh_tott_lht1p0_SC_RI_TCS_SD_C0500'
%  'prs_lh_tott_lht1p0_SC_RI_TCS_SD_C1000'
%  'prs_lh_tott_lht1p0_SC_RI_TCS_SD_C2000'
%  'prs_lh_tott_lht1p0_SC_TR_TCS_SD_C0500'
%  'prs_lh_tott_lht1p0_SC_TR_TCS_SD_C1000'
%  'prs_lh_tott_lht1p0_SC_TR_TCS_SD_C2000'
%  'prs_lh_tott_lht1p0_SC_SS_TCS_SD_C0500'
%  'prs_lh_tott_lht1p0_SC_SS_TCS_SD_C1000'
%  'prs_lh_tott_lht1p0_SC_SS_TCS_SD_C2000'
%  'prs_w_SC_RI_TCS_SD_C0500'
%  'prs_w_SC_RI_TCS_SD_C1000'
%  'prs_w_SC_RI_TCS_SD_C2000'
%  'prs_w_SC_TR_TCS_SD_C0500'
%  'prs_w_SC_TR_TCS_SD_C1000'
%  'prs_w_SC_TR_TCS_SD_C2000'
%  'prs_w_SC_SS_TCS_SD_C0500'
%  'prs_w_SC_SS_TCS_SD_C1000'
%  'prs_w_SC_SS_TCS_SD_C2000'
%  'prs_relvortz_rvz1em4_SC_RI_TCS_SD_C0100'
%  'prs_relvortz_rvz1em4_SC_RI_TCS_SD_C0500'
%  'prs_relvortz_rvz1em4_SC_RI_TCS_SD_C1000'
%  'prs_relvortz_rvz1em4_SC_RI_TCS_SD_C2000'
%  'prs_relvortz_rvz1em4_SC_TR_TCS_SD_C0100'
%  'prs_relvortz_rvz1em4_SC_TR_TCS_SD_C0500'
%  'prs_relvortz_rvz1em4_SC_TR_TCS_SD_C1000'
%  'prs_relvortz_rvz1em4_SC_TR_TCS_SD_C2000'
%  'prs_relvortz_rvz1em4_SC_SS_TCS_SD_C0100'
%  'prs_relvortz_rvz1em4_SC_SS_TCS_SD_C0500'
%  'prs_relvortz_rvz1em4_SC_SS_TCS_SD_C1000'
%  'prs_relvortz_rvz1em4_SC_SS_TCS_SD_C2000'
%  'prs_relvortz_rvz1em4_diff_SC_RI_TCS_SD_C0500'
%  'prs_relvortz_rvz1em4_diff_SC_RI_TCS_SD_C1000'
%  'prs_relvortz_rvz1em4_diff_SC_RI_TCS_SD_C2000'
%  'prs_relvortz_rvz1em4_diff_SC_TR_TCS_SD_C0500'
%  'prs_relvortz_rvz1em4_diff_SC_TR_TCS_SD_C1000'
%  'prs_relvortz_rvz1em4_diff_SC_TR_TCS_SD_C2000'
%  'prs_relvortz_rvz1em4_diff_SC_SS_TCS_SD_C0500'
%  'prs_relvortz_rvz1em4_diff_SC_SS_TCS_SD_C1000'
%  'prs_relvortz_rvz1em4_diff_SC_SS_TCS_SD_C2000'
%  'prs_vaptott_vt0p5_RB_RI_TCS_SD_C0500'
%  'prs_vaptott_vt0p5_RB_RI_TCS_SD_C1000'
%  'prs_vaptott_vt0p5_RB_RI_TCS_SD_C2000'
%  'prs_vaptott_vt0p5_RB_TR_TCS_SD_C0500'
%  'prs_vaptott_vt0p5_RB_TR_TCS_SD_C1000'
%  'prs_vaptott_vt0p5_RB_TR_TCS_SD_C2000'
%  'prs_vaptott_vt0p5_RB_SS_TCS_SD_C0500'
%  'prs_vaptott_vt0p5_RB_SS_TCS_SD_C1000'
%  'prs_vaptott_vt0p5_RB_SS_TCS_SD_C2000'
%  'prs_lh_tott_lht1p0_RB_RI_TCS_SD_C0500'
%  'prs_lh_tott_lht1p0_RB_RI_TCS_SD_C1000'
%  'prs_lh_tott_lht1p0_RB_RI_TCS_SD_C2000'
%  'prs_lh_tott_lht1p0_RB_TR_TCS_SD_C0500'
%  'prs_lh_tott_lht1p0_RB_TR_TCS_SD_C1000'
%  'prs_lh_tott_lht1p0_RB_TR_TCS_SD_C2000'
%  'prs_lh_tott_lht1p0_RB_SS_TCS_SD_C0500'
%  'prs_lh_tott_lht1p0_RB_SS_TCS_SD_C1000'
%  'prs_lh_tott_lht1p0_RB_SS_TCS_SD_C2000'
%  'prs_w_RB_RI_TCS_SD_C0500'
%  'prs_w_RB_RI_TCS_SD_C1000'
%  'prs_w_RB_RI_TCS_SD_C2000'
%  'prs_w_RB_TR_TCS_SD_C0500'
%  'prs_w_RB_TR_TCS_SD_C1000'
%  'prs_w_RB_TR_TCS_SD_C2000'
%  'prs_w_RB_SS_TCS_SD_C0500'
%  'prs_w_RB_SS_TCS_SD_C1000'
%  'prs_w_RB_SS_TCS_SD_C2000'
%  'prs_relhum_RB_SS_TCS_SD_C0500'
%  'prs_relhum_RB_SS_TCS_SD_C1000'
%  'prs_relhum_RB_SS_TCS_SD_C2000'
%  'prs_theta_e_RB_SS_TCS_SD_C0500'
%  'prs_theta_e_RB_SS_TCS_SD_C1000'
%  'prs_theta_e_RB_SS_TCS_SD_C2000'
%  'prs_speed_r_RB_SS_TCS_SD_C0500'
%  'prs_speed_r_RB_SS_TCS_SD_C1000'
%  'prs_speed_r_RB_SS_TCS_SD_C2000'
%  'prs_speed_r_abs_RB_SS_TCS_SD_C0100'
%  'prs_speed_r_abs_RB_SS_TCS_SD_C0500'
%  'prs_speed_r_abs_RB_SS_TCS_SD_C1000'
%  'prs_speed_r_abs_RB_SS_TCS_SD_C2000'
%  };

FigList = {
  'prof_lh_tott_lht1p0_SC_SS'
  'prof_lh_tott_lht1p0_RB_SS'
  'prof_w_SC_SS'
  'prof_w_RB_SS'
  'prof_relhum_SC_SS'
  'prof_relhum_RB_SS'
  'prof_theta_e_SC_SS'
  'prof_theta_e_RB_SS'
  'prof_speed_r_SC_SS'
  'prof_speed_r_RB_SS'
  };


%FigList = {
%  'TcAerosolIngest'
%  };

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
