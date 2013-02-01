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
%            'KeVt_ext'
%            'TS_rmw',
%            'TS_Vint_Vapt_EW'
%            'TS_Vint_Vapt_CO'
%            'TS_Vint_Vapt_RB'
%            'TS_Vint_Lht_EW'
%            'TS_Vint_Lht_CO'
%            'TS_Vint_Lht_RB' };

% FigList = { 'EOF1_theta_e_TCS_SD_C0500'
%             'EOF1_theta_e_TCS_SD_C1000'
%             'EOF1_theta_e_TCS_SD_C2000'
%             'PC1_theta_e_TCS_SD_C0500'
%             'PC1_theta_e_TCS_SD_C1000'
%             'PC1_theta_e_TCS_SD_C2000' };

% FigList = { 'pts_cloud_cm3_twp4_EW_LC_TCS_SD_C0100'
%             'pts_cloud_cm3_twp4_EW_LC_TCS_SD_C0500'
%             'pts_cloud_cm3_twp4_EW_LC_TCS_SD_C1000'
%             'pts_cloud_cm3_twp4_EW_LC_TCS_SD_C2000'
%             'pts_cloud_d_twp4_EW_LC_TCS_SD_C0100'
%             'pts_cloud_d_twp4_EW_LC_TCS_SD_C0500'
%             'pts_cloud_d_twp4_EW_LC_TCS_SD_C1000'
%             'pts_cloud_d_twp4_EW_LC_TCS_SD_C2000'
%             'pts_lh_tott_twp4_lht1p0_EW_LC_TCS_SD_C0100'
%             'pts_lh_tott_twp4_lht1p0_EW_LC_TCS_SD_C0500'
%             'pts_lh_tott_twp4_lht1p0_EW_LC_TCS_SD_C1000'
%             'pts_lh_tott_twp4_lht1p0_EW_LC_TCS_SD_C2000'
%             
%             'pts_cloud_cm3_twp4_CO_LC_TCS_SD_C0100'
%             'pts_cloud_cm3_twp4_CO_LC_TCS_SD_C0500'
%             'pts_cloud_cm3_twp4_CO_LC_TCS_SD_C1000'
%             'pts_cloud_cm3_twp4_CO_LC_TCS_SD_C2000'
%             'pts_cloud_d_twp4_CO_LC_TCS_SD_C0100'
%             'pts_cloud_d_twp4_CO_LC_TCS_SD_C0500'
%             'pts_cloud_d_twp4_CO_LC_TCS_SD_C1000'
%             'pts_cloud_d_twp4_CO_LC_TCS_SD_C2000'
%             'pts_lh_tott_twp4_lht1p0_CO_LC_TCS_SD_C0100'
%             'pts_lh_tott_twp4_lht1p0_CO_LC_TCS_SD_C0500'
%             'pts_lh_tott_twp4_lht1p0_CO_LC_TCS_SD_C1000'
%             'pts_lh_tott_twp4_lht1p0_CO_LC_TCS_SD_C2000'
%             
%             'pts_cloud_cm3_twp4_RB_LC_TCS_SD_C0100'
%             'pts_cloud_cm3_twp4_RB_LC_TCS_SD_C0500'
%             'pts_cloud_cm3_twp4_RB_LC_TCS_SD_C1000'
%             'pts_cloud_cm3_twp4_RB_LC_TCS_SD_C2000'
%             'pts_cloud_d_twp4_RB_LC_TCS_SD_C0100'
%             'pts_cloud_d_twp4_RB_LC_TCS_SD_C0500'
%             'pts_cloud_d_twp4_RB_LC_TCS_SD_C1000'
%             'pts_cloud_d_twp4_RB_LC_TCS_SD_C2000'
%             'pts_lh_tott_twp4_lht1p0_RB_LC_TCS_SD_C0100'
%             'pts_lh_tott_twp4_lht1p0_RB_LC_TCS_SD_C0500'
%             'pts_lh_tott_twp4_lht1p0_RB_LC_TCS_SD_C1000'
%             'pts_lh_tott_twp4_lht1p0_RB_LC_TCS_SD_C2000' };

%FigList = { 'prs_vaptott_twp4_vt0p5_EW_RI_TCS_SD_C0500'
%            'prs_vaptott_twp4_vt0p5_EW_RI_TCS_SD_C1000'
%            'prs_vaptott_twp4_vt0p5_EW_RI_TCS_SD_C2000'
%            'prs_vaptott_twp4_vt0p5_EW_TR_TCS_SD_C0500'
%            'prs_vaptott_twp4_vt0p5_EW_TR_TCS_SD_C1000'
%            'prs_vaptott_twp4_vt0p5_EW_TR_TCS_SD_C2000'
%            'prs_vaptott_twp4_vt0p5_EW_SS_TCS_SD_C0500'
%            'prs_vaptott_twp4_vt0p5_EW_SS_TCS_SD_C1000'
%            'prs_vaptott_twp4_vt0p5_EW_SS_TCS_SD_C2000'
%            'prs_lh_tott_twp4_lht1p0_EW_RI_TCS_SD_C0500'
%            'prs_lh_tott_twp4_lht1p0_EW_RI_TCS_SD_C1000'
%            'prs_lh_tott_twp4_lht1p0_EW_RI_TCS_SD_C2000'
%            'prs_lh_tott_twp4_lht1p0_EW_TR_TCS_SD_C0500'
%            'prs_lh_tott_twp4_lht1p0_EW_TR_TCS_SD_C1000'
%            'prs_lh_tott_twp4_lht1p0_EW_TR_TCS_SD_C2000'
%            'prs_lh_tott_twp4_lht1p0_EW_SS_TCS_SD_C0500'
%            'prs_lh_tott_twp4_lht1p0_EW_SS_TCS_SD_C1000'
%            'prs_lh_tott_twp4_lht1p0_EW_SS_TCS_SD_C2000'
%            'prs_w_twp4_EW_RI_TCS_SD_C0500'
%            'prs_w_twp4_EW_RI_TCS_SD_C1000'
%            'prs_w_twp4_EW_RI_TCS_SD_C2000'
%            'prs_w_twp4_EW_TR_TCS_SD_C0500'
%            'prs_w_twp4_EW_TR_TCS_SD_C1000'
%            'prs_w_twp4_EW_TR_TCS_SD_C2000'
%            'prs_w_twp4_EW_SS_TCS_SD_C0500'
%            'prs_w_twp4_EW_SS_TCS_SD_C1000'
%            'prs_w_twp4_EW_SS_TCS_SD_C2000'
%            'prs_relvortz_rvz1em4_EW_RI_TCS_SD_C0100'
%            'prs_relvortz_rvz1em4_EW_RI_TCS_SD_C0500'
%            'prs_relvortz_rvz1em4_EW_RI_TCS_SD_C1000'
%            'prs_relvortz_rvz1em4_EW_RI_TCS_SD_C2000'
%            'prs_relvortz_rvz1em4_EW_TR_TCS_SD_C0100'
%            'prs_relvortz_rvz1em4_EW_TR_TCS_SD_C0500'
%            'prs_relvortz_rvz1em4_EW_TR_TCS_SD_C1000'
%            'prs_relvortz_rvz1em4_EW_TR_TCS_SD_C2000'
%            'prs_relvortz_rvz1em4_EW_SS_TCS_SD_C0100'
%            'prs_relvortz_rvz1em4_EW_SS_TCS_SD_C0500'
%            'prs_relvortz_rvz1em4_EW_SS_TCS_SD_C1000'
%            'prs_relvortz_rvz1em4_EW_SS_TCS_SD_C2000'
%            'prs_relvortz_rvz1em4_diff_EW_RI_TCS_SD_C0500'
%            'prs_relvortz_rvz1em4_diff_EW_RI_TCS_SD_C1000'
%            'prs_relvortz_rvz1em4_diff_EW_RI_TCS_SD_C2000'
%            'prs_relvortz_rvz1em4_diff_EW_TR_TCS_SD_C0500'
%            'prs_relvortz_rvz1em4_diff_EW_TR_TCS_SD_C1000'
%            'prs_relvortz_rvz1em4_diff_EW_TR_TCS_SD_C2000'
%            'prs_relvortz_rvz1em4_diff_EW_SS_TCS_SD_C0500'
%            'prs_relvortz_rvz1em4_diff_EW_SS_TCS_SD_C1000'
%            'prs_relvortz_rvz1em4_diff_EW_SS_TCS_SD_C2000' };

%FigList = { 
%  'prs_vaptott_twp4_vt0p5_CO_RI_TCS_SD_C0500'
%  'prs_vaptott_twp4_vt0p5_CO_RI_TCS_SD_C1000'
%  'prs_vaptott_twp4_vt0p5_CO_RI_TCS_SD_C2000'
%  'prs_vaptott_twp4_vt0p5_CO_TR_TCS_SD_C0500'
%  'prs_vaptott_twp4_vt0p5_CO_TR_TCS_SD_C1000'
%  'prs_vaptott_twp4_vt0p5_CO_TR_TCS_SD_C2000'
%  'prs_vaptott_twp4_vt0p5_CO_SS_TCS_SD_C0500'
%  'prs_vaptott_twp4_vt0p5_CO_SS_TCS_SD_C1000'
%  'prs_vaptott_twp4_vt0p5_CO_SS_TCS_SD_C2000'
%  'prs_lh_tott_twp4_lht1p0_CO_RI_TCS_SD_C0500'
%  'prs_lh_tott_twp4_lht1p0_CO_RI_TCS_SD_C1000'
%  'prs_lh_tott_twp4_lht1p0_CO_RI_TCS_SD_C2000'
%  'prs_lh_tott_twp4_lht1p0_CO_TR_TCS_SD_C0500'
%  'prs_lh_tott_twp4_lht1p0_CO_TR_TCS_SD_C1000'
%  'prs_lh_tott_twp4_lht1p0_CO_TR_TCS_SD_C2000'
%  'prs_lh_tott_twp4_lht1p0_CO_SS_TCS_SD_C0500'
%  'prs_lh_tott_twp4_lht1p0_CO_SS_TCS_SD_C1000'
%  'prs_lh_tott_twp4_lht1p0_CO_SS_TCS_SD_C2000'
%  'prs_w_twp4_CO_RI_TCS_SD_C0500'
%  'prs_w_twp4_CO_RI_TCS_SD_C1000'
%  'prs_w_twp4_CO_RI_TCS_SD_C2000'
%  'prs_w_twp4_CO_TR_TCS_SD_C0500'
%  'prs_w_twp4_CO_TR_TCS_SD_C1000'
%  'prs_w_twp4_CO_TR_TCS_SD_C2000'
%  'prs_w_twp4_CO_SS_TCS_SD_C0500'
%  'prs_w_twp4_CO_SS_TCS_SD_C1000'
%  'prs_w_twp4_CO_SS_TCS_SD_C2000'
%  };

%FigList = {
%  'prof_lh_tott_twp4_lht1p0_EW_SS'
%  'prof_lh_tott_twp4_lht1p0_CO_SS'
%  'prof_lh_tott_twp4_lht1p0_RB_SS'
%  'prof_w_twp4_EW_SS'
%  'prof_w_twp4_CO_SS'
%  'prof_w_twp4_RB_SS'
%  'prof_relhum_EW_SS'
%  'prof_relhum_CO_SS'
%  'prof_relhum_RB_SS'
%  'prof_theta_e_EW_SS'
%  'prof_theta_e_CO_SS'
%  'prof_theta_e_RB_SS'
%  'prs_relhum_CO_SS_TCS_SD_C0500'
%  'prs_relhum_CO_SS_TCS_SD_C1000'
%  'prs_relhum_CO_SS_TCS_SD_C2000'
%  'prs_theta_e_CO_SS_TCS_SD_C0500'
%  'prs_theta_e_CO_SS_TCS_SD_C1000'
%  'prs_theta_e_CO_SS_TCS_SD_C2000'
%  'vr_diff_TCS_SD_C0500'
%  'vr_diff_TCS_SD_C1000'
%  'vr_diff_TCS_SD_C2000'
%  };

FigList = {
  'TcAerosolIngest'
  };

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
