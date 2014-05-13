function [ ] = GenBarGraphsCtype(ConfigFile)
% GenBarGraphsCtype generate bar plots using cloud type filtered data

  % processing bgraph_<var>.h5 files
  % Two main variables: Averages and Npoints which are organized
  % as: (v,s,c,g) where
  %   v -> Variable name (cot_TSTART, eg)
  %   s -> SST   (ascending order)
  %   c -> CCN   (ascending order)
  %   g -> GCCN  (ascending order)

  % Read the config file to get the structure of how the data is laid out in
  % the file system.
  [ Config ] = ReadConfig(ConfigFile);
    
  Ddir = Config.DiagDir;
  Pdir = Config.PlotDir;

  % Each entry in PlotDefs defines how to build a single plot. The format is:
  %  {
  %  Name
  %  Input file
  %  Title  --> { PanelMarker Label }
  %  Xlabel
  %  Ylabel
  %  LegSpec --> { LegText LegLoc }
  %  Bcolors
  %     entries are names that str2rgb recognizes
  %  Pstyle   ('grouped', 'stacked', etc.)
  %  VarSpec --> { InVarName InVarScale Select Transpose Vname Vmin Vmax }
  %    InVarName is name of dataset inside the HDF5 file
  %    InVarScale is a scaling factor to be applied to InVarName data
  %    Select has to be such that the result will reduce to a 2D array
  %       use vectors (okay to have 1 element vector) for selection
  %       eg: { [1] [1:3] [ 2 4 6 ] [3] } says to select
  %           Averages(1,1:3,[2 4 6],3) reducing it down to 2 dimensions
  %    Then Transpose says whether or not to transpose the array
  %    after it is selected down to two dimensions
  %       1 --> transpose
  %       0 --> do not transpose
  %    Vmin and Vmax are for scaling the y-axis
  %  Output file
  %
  % Need to keep this in sync with GenBarGraphFilesCtype.m
  %   Var order (selection for the first dimension, v, of main variables)
  %     1  _TSTART
  %     2  _TMID
  %     3  _TEND
  %     4  _TALL
  %
  %     5  _strnp_TSTART
  %     6  _strnp_TMID
  %     7  _strnp_TEND
  %     8  _strnp_TALL
  %
  %     9  _strat_TSTART
  %    10  _strat_TMID
  %    11  _strat_TEND
  %    12  _strat_TALL
  %
  %    13  _cumul_TSTART
  %    14  _cumul_TMID
  %    15  _cumul_TEND
  %    16  _cumul_TALL
  %
  %    17  _all_cld_TSTART
  %    18  _all_cld_TMID
  %    19  _all_cld_TEND
  %    20  _all_cld_TALL
  %
  %    21  _stall_TSTART
  %    22  _stall_TMID
  %    23  _stall_TEND
  %    24  _stall_TALL
  %
  % Cloud fraction is an exception to this rule, where there is no _all_cld. Instead
  % there is an stmix and scmix case:
  %
  %    ...
  %
  %    17 _stmix_TSTART
  %    18 _stmix_TMID
  %    19 _stmix_TEND
  %    20 _stmix_TALL
  %
  %    21 _scmix_TSTART
  %    22 _scmix_TMID
  %    23 _scmix_TEND
  %    24 _scmix_TALL
  %
  % Added LCL, these are in the bgraph_lcl fiels which have variables ordered as:
  % 
  %     1 lcl_TSTART
  %     2 lcl_TMID
  %     3 lcl_TEND
  %     4 lcl_TALL
  %
  %     5 lcl_stall_TSTART
  %     6 lcl_stall_TMID
  %     7 lcl_stall_TEND
  %     8 lcl_stall_TALL
  %
  % Added entrainment velocities, these are in the bgraph_we files which have
  % variables ordered as:
  %
  %     1 ThetaWe_turb_all_TALL
  %     2 ThetaV_We_turb_all_TALL
  %     3 VaporWe_turb_all_TALL
  %
  %     4 ThetaWe_gm_all_TALL
  %     5 ThetaV_We_gm_all_TALL
  %     6 VaporWe_gm_all_TALL
  %
  %     7 ThetaWe_gm_stall_TALL
  %     8 ThetaV_We_gm_stall_TALL
  %     9 VaporWe_gm_stall_TALL
  %

  PlotDefs = {
       %%%%%%%%%%%% CLOUD OPTICAL THICKNESS %%%%%%%%%%%%%%%%%%%
       % COT averages, domain, all time points, CCN only
       {
       'COT Avg, DOMAIN TALL CCN only, grouped by SST'
       'bgraph_cot.h5'
       { 'a' 'Domain' }
       'N_a (# cm^-^3)'
       '\tau_c'
       { 'blue' 'cyan' 'magenta' }
       { { 'S293', 'S298', 'S303' } 'NorthWest' }
       'grouped'
       { 'Averages' 1 { [4] [1:3] [1:6] [1] } 1 'CCN' 0 4.5 }
       'bars_avg_cot_TALL_CO.jpg'
       }

       % COT averages, all clouds, all time points, CCN only
       {
       'COT Avg, ALL_CLD TALL CCN only, grouped by SST'
       'bgraph_cot.h5'
       { 'b' 'All Clouds' }
       'N_a (# cm^-^3)'
       '\tau_c'
       { 'blue' 'cyan' 'magenta' }
       { { 'S293', 'S298', 'S303' } 'NorthWest' }
       'grouped'
       { 'Averages' 1 { [20] [1:3] [1:6] [1] } 1 'CCN' 0 6 }
       'bars_avg_cot_all_cld_TALL_CO.jpg'
       }

       % COT averages, stall, all time points, CCN only
       {
       'COT Avg, ST TALL CCN only, grouped by SST'
       'bgraph_cot.h5'
       { 'a' 'ST' }
       'N_a (# cm^-^3)'
       '\tau_c'
       { 'blue' 'cyan' 'magenta' }
       { { 'S293', 'S298', 'S303' } 'NorthWest' }
       'grouped'
       { 'Averages' 1 { [24] [1:3] [1:6] [1] } 1 'CCN' 0 4.5 }
       'bars_avg_cot_stall_TALL_CO.jpg'
       }

       %%%%%%%%%%%% LIQUID WATER PATH %%%%%%%%%%%%%%%%%%%
       % LWP averages, domain, all time points, CCN only
       {
       'LWP Avg, Domain, TALL CCN only, grouped by SST'
       'bgraph_lwp.h5'
       { 'a' 'Domain' }
       'N_a (# cm^-^3)'
       'LWP (mm)'
       { 'blue' 'cyan' 'magenta' }
       { { 'S293', 'S298', 'S303' } 'NorthWest' }
       'grouped'
       { 'Averages' 1 { [4] [1:3] [1:6] [1] } 1 'CCN' 0 0.25 }
       'bars_avg_lwp_TALL_CO.jpg'
       }

       % LWP averages, all clouds, time points, CCN only
       {
       'LWP Avg, ALL_CLD, TALL CCN only, grouped by SST'
       'bgraph_lwp.h5'
       { 'a' 'All Clouds' }
       'N_a (# cm^-^3)'
       'LWP (mm)'
       { 'blue' 'cyan' 'magenta' }
       { { 'S293', 'S298', 'S303' } 'NorthWest' }
       'grouped'
       { 'Averages' 1 { [20] [1:3] [1:6] [1] } 1 'CCN' 0 0.25 }
       'bars_avg_lwp_all_cld_TALL_CO.jpg'
       }

       % LWP averages, stall, time points, CCN only
       {
       'LWP Avg, ST, TALL CCN only, grouped by SST'
       'bgraph_lwp.h5'
       { 'a' 'ST' }
       'N_a (# cm^-^3)'
       'LWP (mm)'
       { 'blue' 'cyan' 'magenta' }
       { { 'S293', 'S298', 'S303' } 'NorthWest' }
       'grouped'
       { 'Averages' 1 { [24] [1:3] [1:6] [1] } 1 'CCN' 0 0.25 }
       'bars_avg_lwp_stall_TALL_CO.jpg'
       }

       %%%%%%%%%%%%%% CLOUD DEPTH %%%%%%%%%%%%%%%%%%%%%%%%
       % Cloud depth averages, domain, all time points, CCN only
       {
       'Cloud Depth Avg, Domain, TALL CCN only, grouped by SST'
       'bgraph_cdepth.h5'
       { 'a' 'Domain' }
       'N_a (# cm^-^3)'
       'Cloud Depth (m)'
       { 'blue' 'cyan' 'magenta' }
       { { 'S293', 'S298', 'S303' } 'NorthWest' }
       'grouped'
       { 'Averages' 1 { [4] [1:3] [1:6] [1] } 1 'CCN' 0 900 }
       'bars_avg_cdepth_TALL_CO.jpg'
       }

       % Cloud depth averages, all clouds, all time points, CCN only
       {
       'Cloud Depth Avg, ALL_CLD, TALL CCN only, grouped by SST'
       'bgraph_cdepth.h5'
       { 'a' 'All Clouds' }
       'N_a (# cm^-^3)'
       'Cloud Depth (m)'
       { 'blue' 'cyan' 'magenta' }
       { { 'S293', 'S298', 'S303' } 'NorthWest' }
       'grouped'
       { 'Averages' 1 { [20] [1:3] [1:6] [1] } 1 'CCN' 0 900 }
       'bars_avg_cdepth_all_cld_TALL_CO.jpg'
       }

       % Cloud depth averages, stall, all time points, CCN only
       {
       'Cloud Depth Avg, stall, TALL CCN only, grouped by SST'
       'bgraph_cdepth.h5'
       { 'a' 'ST' }
       'N_a (# cm^-^3)'
       'Cloud Depth (m)'
       { 'blue' 'cyan' 'magenta' }
       { { 'S293', 'S298', 'S303' } 'NorthWest' }
       'grouped'
       { 'Averages' 1 { [24] [1:3] [1:6] [1] } 1 'CCN' 0 900 }
       'bars_avg_cdepth_stall_TALL_CO.jpg'
       }

       %%%%%%%%%%%%%%%% PRECIP RATE %%%%%%%%%%%%%%%%%%%
       % PR averages, domain, all time points, CCN only
       {
       'PR Avg, Domain, TALL CCN only, grouped by SST'
       'bgraph_pcprr.h5'
       { 'a' 'Domain' }
       'N_a (# cm^-^3)'
       'Precip Rate (mm h^-^1)'
       { 'blue' 'cyan' 'magenta' }
       { { 'S293', 'S298', 'S303' } 'NorthWest' }
       'grouped'
       { 'Averages' 1 { [4] [1:3] [1:6] [1] } 1 'CCN' 0 0.12 }
       'bars_avg_pcprr_TALL_CO.jpg'
       }

       % PR averages, all clouds, all time points, CCN only
       {
       'PR Avg, ALL_CLD, TALL CCN only, grouped by SST'
       'bgraph_pcprr.h5'
       { 'a' 'All Clouds' }
       'N_a (# cm^-^3)'
       'Precip Rate (mm h^-^1)'
       { 'blue' 'cyan' 'magenta' }
       { { 'S293', 'S298', 'S303' } 'NorthWest' }
       'grouped'
       { 'Averages' 1 { [20] [1:3] [1:6] [1] } 1 'CCN' 0 0.3 }
       'bars_avg_pcprr_all_cld_TALL_CO.jpg'
       }

       % PR averages, stall, all time points, CCN only
       {
       'PR Avg, Stall, TALL CCN only, grouped by SST'
       'bgraph_pcprr.h5'
       { 'a' 'ST' }
       'N_a (# cm^-^3)'
       'Precip Rate (mm h^-^1)'
       { 'blue' 'cyan' 'magenta' }
       { { 'S293', 'S298', 'S303' } 'NorthWest' }
       'grouped'
       { 'Averages' 1 { [24] [1:3] [1:6] [1] } 1 'CCN' 0 0.07 }
       'bars_avg_pcprr_stall_TALL_CO.jpg'
       }


%       % LCL averages, DOMAIN all time points, CCN only
%       {
%       'LCL Avg TALL CCN only, grouped by SST'
%       'bgraph_lcl.h5'
%       { 'a' 'Domain' }
%       'N_a (# cm^-^3)'
%       'LCL (m)'
%       { 'blue' 'cyan' 'magenta' }
%       { { 'S293', 'S298', 'S303' } 'NorthWest' }
%       'grouped'
%       { 'Averages' 1 { [4] [1:3] [1:6] [1] } 1 'CCN' 0 1000 }
%       'bars_avg_lcl_TALL_CO.jpg'
%       }
%
%       % Theta We averages, all time points, CCN only
%       {
%       'We theta Avg, DOMAIN TALL CCN only, grouped by SST'
%       'bgraph_we.h5'
%       { 'a' 'Domain' }
%       'N_a (# cm^-^3)'
%       '\theta W_e (m s^-^1)'
%       { 'blue' 'cyan' 'magenta' }
%       { { 'S293', 'S298', 'S303' } 'NorthWest' }
%       'grouped'
%       { 'Averages' 1 { [4] [1:3] [1:6] [1] } 1 'CCN' 0 0.03 }
%       'bars_avg_theta_we_TALL_CO.jpg'
%       }
%
%       % ThetaV We averages, all time points, CCN only
%       {
%       'We theta_v Avg, DOMAIN TALL CCN only, grouped by SST'
%       'bgraph_we.h5'
%       { 'a' 'Domain' }
%       'N_a (# cm^-^3)'
%       '\theta_v W_e (m s^-^1)'
%       { 'blue' 'cyan' 'magenta' }
%       { { 'S293', 'S298', 'S303' } 'NorthWest' }
%       'grouped'
%       { 'Averages' 1 { [5] [1:3] [1:6] [1] } 1 'CCN' 0 0.06 }
%       'bars_avg_theta_v_we_TALL_CO.jpg'
%       }
%
%       % Vapor We averages, all time points, CCN only
%       {
%       'We vapor Avg, DOMAIN TALL CCN only, grouped by SST'
%       'bgraph_we.h5'
%       { 'a' 'Domain' }
%       'N_a (# cm^-^3)'
%       'q_v W_e (m s^-^1)'
%       { 'blue' 'cyan' 'magenta' }
%       { { 'S293', 'S298', 'S303' } 'NorthWest' }
%       'grouped'
%       { 'Averages' 1 { [6] [1:3] [1:6] [1] } 1 'CCN' -0.1 0.05 }
%       'bars_avg_vapor_we_TALL_CO.jpg'
%       }
%
%
%       % LCL averages, ST all time points, CCN only
%       {
%       'LCL Avg TALL CCN only, grouped by SST'
%       'bgraph_lcl.h5'
%       { 'a' 'ST' }
%       'N_a (# cm^-^3)'
%       'LCL (m)'
%       { 'blue' 'cyan' 'magenta' }
%       { { 'S293', 'S298', 'S303' } 'NorthWest' }
%       'grouped'
%       { 'Averages' 1 { [8] [1:3] [1:6] [1] } 1 'CCN' 0 1000 }
%       'bars_avg_lcl_stall_TALL_CO.jpg'
%       }
%
%       % LCL averages, ST all time points, CCN only, S293
%       {
%       'LCL Avg TALL CCN only, S293'
%       'bgraph_lcl.h5'
%       { 'a' 'ST, S293' }
%       'N_a (# cm^-^3)'
%       'LCL (m)'
%       { 'cyan' }
%       { { 'NoLegend' } 'NorthWest' }
%       'grouped'
%       { 'Averages' 1 { [8] [1] [1:6] [1] } 1 'CCN' 0 400 }
%       'bars_avg_lcl_stall_TALL_CO_S293.jpg'
%       }
%
%       % LCL averages, ST all time points, CCN only, S298
%       {
%       'LCL Avg TALL CCN only, S298'
%       'bgraph_lcl.h5'
%       { 'a' 'ST, S298' }
%       'N_a (# cm^-^3)'
%       'LCL (m)'
%       { 'cyan' }
%       { { 'NoLegend' } 'NorthWest' }
%       'grouped'
%       { 'Averages' 1 { [8] [2] [1:6] [1] } 1 'CCN' 0 600 }
%       'bars_avg_lcl_stall_TALL_CO_S298.jpg'
%       }
%
%       % LCL averages, ST all time points, CCN only, S303
%       {
%       'LCL Avg TALL CCN only, S303'
%       'bgraph_lcl.h5'
%       { 'a' 'ST, S303' }
%       'N_a (# cm^-^3)'
%       'LCL (m)'
%       { 'cyan' }
%       { { 'NoLegend' } 'NorthWest' }
%       'grouped'
%       { 'Averages' 1 { [8] [3] [1:6] [1] } 1 'CCN' 0 900 }
%       'bars_avg_lcl_stall_TALL_CO_S303.jpg'
%       }
%
%       % Theta We averages, all time points, CCN only
%       {
%       'We theta Avg, ST TALL CCN only, grouped by SST'
%       'bgraph_we.h5'
%       { 'a' 'ST' }
%       'N_a (# cm^-^3)'
%       '\theta W_e (m s^-^1)'
%       { 'blue' 'cyan' 'magenta' }
%       { { 'S293', 'S298', 'S303' } 'NorthWest' }
%       'grouped'
%       { 'Averages' 1 { [7] [1:3] [1:6] [1] } 1 'CCN' 0 0.012 }
%       'bars_avg_theta_we_stall_TALL_CO.jpg'
%       }
%
%       % ThetaV We averages, all time points, CCN only
%       {
%       'We theta_v Avg, ST TALL CCN only, grouped by SST'
%       'bgraph_we.h5'
%       { 'a' 'ST' }
%       'N_a (# cm^-^3)'
%       '\theta_v W_e (cm s^-^1)'
%       { 'blue' 'cyan' 'magenta' }
%       { { 'S293', 'S298', 'S303' } 'NorthWest' }
%       'grouped'
%       { 'Averages' 100 { [8] [1:3] [1:6] [1] } 1 'CCN' 0 1.2 }
%       'bars_avg_theta_v_we_stall_TALL_CO.jpg'
%       }
%
%       % ThetaV We averages, all time points, CCN only, S293
%       {
%       'We theta_v Avg, ST TALL CCN only, S293'
%       'bgraph_we.h5'
%       { 'f' 'ST, S293' }
%       'N_a (# cm^-^3)'
%       '\theta_v W_e (cm s^-^1)'
%       { 'cyan' }
%       { { 'NoLegend' } 'NorthWest' }
%       'grouped'
%       { 'Averages' 100 { [8] [1] [1:6] [1] } 1 'CCN' 0 0.15 }
%       'bars_avg_theta_v_we_stall_TALL_CO_S293.jpg'
%       }
%
%       % ThetaV We averages, all time points, CCN only, S298
%       {
%       'We theta_v Avg, ST TALL CCN only, S298'
%       'bgraph_we.h5'
%       { 'f' 'ST, S298' }
%       'N_a (# cm^-^3)'
%       '\theta_v W_e (cm s^-^1)'
%       { 'cyan' }
%       { { 'NoLegend' } 'NorthWest' }
%       'grouped'
%       { 'Averages' 100 { [8] [2] [1:6] [1] } 1 'CCN' 0 0.6 }
%       'bars_avg_theta_v_we_stall_TALL_CO_S298.jpg'
%       }
%
%       % ThetaV We averages, all time points, CCN only, S303
%       {
%       'We theta_v Avg, ST TALL CCN only, S303'
%       'bgraph_we.h5'
%       { 'f' 'ST, S303' }
%       'N_a (# cm^-^3)'
%       '\theta_v W_e (cm s^-^1)'
%       { 'cyan' }
%       { { 'NoLegend' } 'NorthWest' }
%       'grouped'
%       { 'Averages' 100 { [8] [3] [1:6] [1] } 1 'CCN' 0 1.1 }
%       'bars_avg_theta_v_we_stall_TALL_CO_S303.jpg'
%       }
%
%       % ThetaV We averages, all time points, CCN only, C50
%       {
%       'We theta_v Avg, ST TALL CCN only, C50'
%       'bgraph_we.h5'
%       { 'f' 'ST, C50' }
%       'SST (K)'
%       '\theta_v W_e (cm s^-^1)'
%       { 'cyan' }
%       { { 'NoLegend' } 'NorthWest' }
%       'grouped'
%       { 'Averages' 100 { [8] [1:3] [1] [1] } 1 'SST' 0 1.1 }
%       'bars_avg_theta_v_we_stall_TALL_CO_C50.jpg'
%       }
%
%       % ThetaV We averages, all time points, CCN only, C400
%       {
%       'We theta_v Avg, ST TALL CCN only, C400'
%       'bgraph_we.h5'
%       { 'f' 'ST, C400' }
%       'SST (K)'
%       '\theta_v W_e (cm s^-^1)'
%       { 'cyan' }
%       { { 'NoLegend' } 'NorthWest' }
%       'grouped'
%       { 'Averages' 100 { [8] [1:3] [4] [1] } 1 'SST' 0 0.6 }
%       'bars_avg_theta_v_we_stall_TALL_CO_C400.jpg'
%       }
%
%       % ThetaV We averages, all time points, CCN only, C1600
%       {
%       'We theta_v Avg, ST TALL CCN only, C1600'
%       'bgraph_we.h5'
%       { 'f' 'ST, C1600' }
%       'SST (K)'
%       '\theta_v W_e (cm s^-^1)'
%       { 'cyan' }
%       { { 'NoLegend' } 'NorthWest' }
%       'grouped'
%       { 'Averages' 100 { [8] [1:3] [6] [1] } 1 'SST' 0 0.6 }
%       'bars_avg_theta_v_we_stall_TALL_CO_C1600.jpg'
%       }
%
%       % Vapor We averages, all time points, CCN only
%       {
%       'We vapor Avg, ST TALL CCN only, grouped by SST'
%       'bgraph_we.h5'
%       { 'a' 'ST' }
%       'N_a (# cm^-^3)'
%       'q_v W_e (cm s^-^1)'
%       { 'blue' 'cyan' 'magenta' }
%       { { 'S293', 'S298', 'S303' } 'NorthWest' }
%       'grouped'
%       { 'Averages' 100 { [9] [1:3] [1:6] [1] } 1 'CCN' -0.5 12 }
%       'bars_avg_vapor_we_stall_TALL_CO.jpg'
%       }
%
%       % Vapor We averages, all time points, CCN only, S293
%       {
%       'We vapor Avg, ST TALL CCN only, S293'
%       'bgraph_we.h5'
%       { 'h' 'ST, S293' }
%       'N_a (# cm^-^3)'
%       'q_v W_e (cm s^-^1)'
%       { 'cyan' }
%       { { 'NoLegend' } 'NorthWest' }
%       'grouped'
%       { 'Averages' 100 { [9] [1] [1:6] [1] } 1 'CCN' 0 0.1 }
%       'bars_avg_vapor_we_stall_TALL_CO_S293.jpg'
%       }
%
%       % Vapor We averages, all time points, CCN only, S298
%       {
%       'We vapor Avg, ST TALL CCN only, S298'
%       'bgraph_we.h5'
%       { 'h' 'ST, S298' }
%       'N_a (# cm^-^3)'
%       'q_v W_e (cm s^-^1)'
%       { 'cyan' }
%       { { 'NoLegend' } 'NorthWest' }
%       'grouped'
%       { 'Averages' 100 { [9] [2] [1:6] [1] } 1 'CCN' 0 0.8 }
%       'bars_avg_vapor_we_stall_TALL_CO_S298.jpg'
%       }
%
%       % Vapor We averages, all time points, CCN only, S303
%       {
%       'We vapor Avg, ST TALL CCN only, S303'
%       'bgraph_we.h5'
%       { 'h' 'ST, S303' }
%       'N_a (# cm^-^3)'
%       'q_v W_e (cm s^-^1)'
%       { 'cyan' }
%       { { 'NoLegend' } 'NorthWest' }
%       'grouped'
%       { 'Averages' 100 { [9] [3] [1:6] [1] } 1 'CCN' -1 12 }
%       'bars_avg_vapor_we_stall_TALL_CO_S303.jpg'
%       }
%
%       % Vapor We averages, all time points, CCN only, C50
%       {
%       'We vapor Avg, ST TALL CCN only, C50'
%       'bgraph_we.h5'
%       { 'h' 'ST, C50' }
%       'SST (K)'
%       'q_v W_e (cm s^-^1)'
%       { 'cyan' }
%       { { 'NoLegend' } 'NorthWest' }
%       'grouped'
%       { 'Averages' 100 { [9] [1:3] [1] [1] } 1 'SST' -0.4 0.4 }
%       'bars_avg_vapor_we_stall_TALL_CO_C50.jpg'
%       }
%
%       % Vapor We averages, all time points, CCN only, C400
%       {
%       'We vapor Avg, ST TALL CCN only, C400'
%       'bgraph_we.h5'
%       { 'h' 'ST, C400' }
%       'SST (K)'
%       'q_v W_e (cm s^-^1)'
%       { 'cyan' }
%       { { 'NoLegend' } 'NorthWest' }
%       'grouped'
%       { 'Averages' 100 { [9] [1:3] [4] [1] } 1 'SST' 0 4 }
%       'bars_avg_vapor_we_stall_TALL_CO_C400.jpg'
%       }
%
%       % Vapor We averages, all time points, CCN only, C1600
%       {
%       'We vapor Avg, ST TALL CCN only, C1600'
%       'bgraph_we.h5'
%       { 'h' 'ST, C1600' }
%       'SST (K)'
%       'q_v W_e (cm s^-^1)'
%       { 'cyan' }
%       { { 'NoLegend' } 'NorthWest' }
%       'grouped'
%       { 'Averages' 100 { [9] [1:3] [6] [1] } 1 'SST' 0 4.5 }
%       'bars_avg_vapor_we_stall_TALL_CO_C1600.jpg'
%       }
%
%       % In subsequent code, extra actions are taken according to the first
%       % entry here containing "Cloud Distribution"
%       %    The first three columns are added together to form an ST type
%       %         (SNP + STM + ST)
%       %    The result is subtracted from 100 to form a Clear type
%       % Cloud type, npoints (relative amounts), TALL, CCN only, SST 293
%       {
%       'Cloud Distribution TALL CCN only, 293K'
%       'bgraph_cfrac.h5'
%       { 'a' 'Domain, S293' }
%       'N_a (# cm^-^3)'
%       'Cloud Distribution (%)'
%       { 'cyan' 'blue' 'magenta' 'white' }
%       { { 'ST' 'SC' 'CV' 'Clear' } 'NorthWest' }
%       'stacked'
%       { 'Averages' 100 { [8 20 12 24 16] [1] [1:6] [1] } 1 'CCN' 0 120 }
%       'bars_avg_ctype_TALL_CO_S293.jpg'
%       }
%
%       % Cloud type, npoints (relative amounts), TALL, CCN only, SST 298
%       {
%       'Cloud Distribution TALL CCN only, 298K'
%       'bgraph_cfrac.h5'
%       { 'b' 'Domain, S298' }
%       'N_a (# cm^-^3)'
%       'Cloud Distribution (%)'
%       { 'cyan' 'blue' 'magenta' 'white' }
%       { { 'ST' 'SC' 'CV' 'Clear' } 'NorthWest' }
%       'stacked'
%       { 'Averages' 100 { [8 20 12 24 16] [2] [1:6] [1] } 1 'CCN' 0 120 }
%       'bars_avg_ctype_TALL_CO_S298.jpg'
%       }
%
%       % Cloud type, npoints (relative amounts), TALL, CCN only, SST 303
%       {
%       'Cloud Distribution TALL CCN only, 300K'
%       'bgraph_cfrac.h5'
%       { 'c' 'Domain, S303' }
%       'N_a (# cm^-^3)'
%       'Cloud Distribution (%)'
%       { 'cyan' 'blue' 'magenta' 'white' }
%       { { 'ST' 'SC' 'CV' 'Clear' } 'NorthWest' }
%       'stacked'
%       { 'Averages' 100 { [8 20 12 24 16] [3] [1:6] [1] } 1 'CCN' 0 120 }
%       'bars_avg_ctype_TALL_CO_S303.jpg'
%       }

%%%       % Stratiform (non-precipitating) cloud types, COT averages, all time points, CCN only
%%%       {
%%%       'Stratiform (NP) Avg TALL CCN only, grouped by SST'
%%%       'bgraph_cot.h5'
%%%       { 'a' 'SNP' }
%%%       'N_a (# cm^-^3)'
%%%       '\tau_c'
%%%       { 'blue' 'cyan' 'magenta' }
%%%       { { 'S293', 'S298', 'S303' } 'NorthWest' }
%%%       'grouped'
%%%       { 'Averages' 1 { [8] [1:3] [1:6] [1] } 1 'CCN' 0 1.5 }
%%%       'bars_avg_cot_strnp_TALL_CO.jpg'
%%%       }
%%%
%%%       % Stratiform (non-precipitating) cloud types, PR averages, all time points, CCN only
%%%       {
%%%       'Stratiform (NP) PR Avg TALL CCN only, grouped by SST'
%%%       'bgraph_pcprr.h5'
%%%       { 'a' 'SNP' }
%%%       'N_a (# cm^-^3)'
%%%       'Precip Rate (mm h^-^1)'
%%%       { 'blue' 'cyan' 'magenta' }
%%%       { { 'S293', 'S298', 'S303' } 'NorthWest' }
%%%       'grouped'
%%%       { 'Averages' 1 { [8] [1:3] [1:6] [1] } 1 'CCN' 0 0.001 }
%%%       'bars_avg_pcprr_strnp_TALL_CO.jpg'
%%%       }
%%%
%%%       % Stratiform (non-precipitating) cloud types, LWP averages, all time points, CCN only
%%%       {
%%%       'Stratiform (NP) LWP Avg TALL CCN only, grouped by SST'
%%%       'bgraph_lwp.h5'
%%%       { 'a' 'SNP' }
%%%       'N_a (# cm^-^3)'
%%%       'LWP (mm)'
%%%       { 'blue' 'cyan' 'magenta' }
%%%       { { 'S293', 'S298', 'S303' } 'NorthWest' }
%%%       'grouped'
%%%       { 'Averages' 1 { [8] [1:3] [1:6] [1] } 1 'CCN' 0 0.035 }
%%%       'bars_avg_lwp_strnp_TALL_CO.jpg'
%%%       }
%%%
%%%       % Stratiform (non-precipitating) cloud types, Cloud depth averages, all time points, CCN only
%%%       {
%%%       'Stratiform (NP) Cloud Depth Avg TALL CCN only, grouped by SST'
%%%       'bgraph_cdepth.h5'
%%%       { 'a' 'SNP' }
%%%       'N_a (# cm^-^3)'
%%%       'Cloud Depth (m)'
%%%       { 'blue' 'cyan' 'magenta' }
%%%       { { 'S293', 'S298', 'S303' } 'NorthWest' }
%%%       'grouped'
%%%       { 'Averages' 1 { [8] [1:3] [1:6] [1] } 1 'CCN' 0 180 }
%%%       'bars_avg_cdepth_strnp_TALL_CO.jpg'
%%%       }
%%%
%%%       % Stratiform (precipitating) cloud types, COT averages, all time points, CCN only
%%%       {
%%%       'Stratiform Avg TALL CCN only, grouped by SST'
%%%       'bgraph_cot.h5'
%%%       { 'a' 'ST' }
%%%       'N_a (# cm^-^3)'
%%%       '\tau_c'
%%%       { 'blue' 'cyan' 'magenta' }
%%%       { { 'S293', 'S298', 'S303' } 'NorthWest' }
%%%       'grouped'
%%%       { 'Averages' 1 { [12] [1:3] [1:6] [1] } 1 'CCN' 0 5 }
%%%       'bars_avg_cot_strat_TALL_CO.jpg'
%%%       }
%%%
%%%       % Stratiform (precipitating) cloud types, PR averages, all time points, CCN only
%%%       {
%%%       'Stratiform PR Avg TALL CCN only, grouped by SST'
%%%       'bgraph_pcprr.h5'
%%%       { 'a' 'ST' }
%%%       'N_a (# cm^-^3)'
%%%       'Precip Rate (mm h^-^1)'
%%%       { 'blue' 'cyan' 'magenta' }
%%%       { { 'S293', 'S298', 'S303' } 'NorthWest' }
%%%       'grouped'
%%%       { 'Averages' 1 { [12] [1:3] [1:6] [1] } 1 'CCN' 0 0.15 }
%%%       'bars_avg_pcprr_strat_TALL_CO.jpg'
%%%       }
%%%
%%%       % Stratiform (precipitating) cloud types, LWP averages, all time points, CCN only
%%%       {
%%%       'Stratiform LWP Avg TALL CCN only, grouped by SST'
%%%       'bgraph_lwp.h5'
%%%       { 'a' 'ST' }
%%%       'N_a (# cm^-^3)'
%%%       'LWP (mm)'
%%%       { 'blue' 'cyan' 'magenta' }
%%%       { { 'S293', 'S298', 'S303' } 'NorthWest' }
%%%       'grouped'
%%%       { 'Averages' 1 { [12] [1:3] [1:6] [1] } 1 'CCN' 0 0.5 }
%%%       'bars_avg_lwp_strat_TALL_CO.jpg'
%%%       }
%%%
%%%       % Stratiform (precipitating) cloud types, Cloud depth averages, all time points, CCN only
%%%       {
%%%       'Stratiform Cloud Depth Avg TALL CCN only, grouped by SST'
%%%       'bgraph_cdepth.h5'
%%%       { 'a' 'ST' }
%%%       'N_a (# cm^-^3)'
%%%       'Cloud Depth (m)'
%%%       { 'blue' 'cyan' 'magenta' }
%%%       { { 'S293', 'S298', 'S303' } 'NorthWest' }
%%%       'grouped'
%%%       { 'Averages' 1 { [12] [1:3] [1:6] [1] } 1 'CCN' 0 900 }
%%%       'bars_avg_cdepth_strat_TALL_CO.jpg'
%%%       }
%%%
%%%       % Convective cloud types, COT averages, all time points, CCN only
%%%       {
%%%       'Convective Avg TALL CCN only, grouped by SST'
%%%       'bgraph_cot.h5'
%%%       { 'a' 'CV' }
%%%       'N_a (# cm^-^3)'
%%%       '\tau_c'
%%%       { 'blue' 'cyan' 'magenta' }
%%%       { { 'S293', 'S298', 'S303' } 'NorthWest' }
%%%       'grouped'
%%%       { 'Averages' 1 { [16] [1:3] [1:6] [1] } 1 'CCN' 0 120 }
%%%       'bars_avg_cot_cumul_TALL_CO.jpg'
%%%       }
%%%
%%%       % Convective cloud types, PR averages, all time points, CCN only
%%%       {
%%%       'Convective PR Avg TALL CCN only, grouped by SST'
%%%       'bgraph_pcprr.h5'
%%%       { 'a' 'CV' }
%%%      'N_a (# cm^-^3)'
%%%       'Precip Rate (mm h^-^1)'
%%%       { 'blue' 'cyan' 'magenta' }
%%%       { { 'S293', 'S298', 'S303' } 'NorthWest' }
%%%       'grouped'
%%%       { 'Averages' 1 { [16] [1:3] [1:6] [1] } 1 'CCN' 0 30 }
%%%       'bars_avg_pcprr_cumul_TALL_CO.jpg'
%%%       }
%%%
%%%       % Convective cloud types, LWP averages, all time points, CCN only
%%%       {
%%%       'Convective LWP Avg TALL CCN only, grouped by SST'
%%%       'bgraph_lwp.h5'
%%%       { 'a' 'CV' }
%%%       'N_a (# cm^-^3)'
%%%       'LWP (mm)'
%%%       { 'blue' 'cyan' 'magenta' }
%%%       { { 'S293', 'S298', 'S303' } 'NorthWest' }
%%%       'grouped'
%%%       { 'Averages' 1 { [16] [1:3] [1:6] [1] } 1 'CCN' 0 5 }
%%%       'bars_avg_lwp_cumul_TALL_CO.jpg'
%%%       }
%%%
%%%       % Convective cloud types, Cloud depth averages, all time points, CCN only
%%%       {
%%%       'Convective Cloud Depth Avg TALL CCN only, grouped by SST'
%%%       'bgraph_cdepth.h5'
%%%       { 'a' 'CV' }
%%%       'N_a (# cm^-^3)'
%%%       'Cloud Depth (m)'
%%%       { 'blue' 'cyan' 'magenta' }
%%%       { { 'S293', 'S298', 'S303' } 'NorthWest' }
%%%       'grouped'
%%%       { 'Averages' 1 { [16] [1:3] [1:6] [1] } 1 'CCN' 0 4000 }
%%%       'bars_avg_cdepth_cumul_TALL_CO.jpg'
%%%       }
%%%
%%%       % All cloud types, COT averages, all time points, CCN only
%%%       {
%%%       'All Clouds Avg TALL CCN only, grouped by SST'
%%%       'bgraph_cot.h5'
%%%       { 'a' 'All Clouds' }
%%%       'N_a (# cm^-^3)'
%%%       '\tau_c'
%%%       { 'blue' 'cyan' 'magenta' }
%%%       { { 'S293', 'S298', 'S303' } 'NorthWest' }
%%%       'grouped'
%%%       { 'Averages' 1 { [20] [1:3] [1:6] [1] } 1 'CCN' 0 6 }
%%%       'bars_avg_cot_all_cld_TALL_CO.jpg'
%%%       }
%%%
%%%       % Cloud fraction averages, all time points, CCN only
%%%       {
%%%       'Cloud Fraction Avg TALL CCN only, grouped by SST'
%%%       'bgraph_cfrac.h5'
%%%       { 'a' 'Domain' }
%%%       'N_a (# cm^-^3)'
%%%       'Cloud Fraction'
%%%       { 'blue' 'cyan' 'magenta' }
%%%       { { 'S293', 'S298', 'S303' } 'NorthWest' }
%%%       'grouped'
%%%       { 'Averages' 1 { [4] [1:3] [1:6] [1] } 1 'CCN' 0 1.2 }
%%%       'bars_avg_cfrac_TALL_CO.jpg'
%%%       }
%%%
%%%       % Stratiform (NP) Cloud fraction averages, all time points, CCN only
%%%       {
%%%       'Stratiform (NP) Cloud Fraction Avg TALL CCN only, grouped by SST'
%%%       'bgraph_cfrac.h5'
%%%       { 'a' 'SNP' }
%%%       'N_a (# cm^-^3)'
%%%       'Cloud Fraction'
%%%       { 'blue' 'cyan' 'magenta' }
%%%       { { 'S293', 'S298', 'S303' } 'NorthWest' }
%%%       'grouped'
%%%       { 'Averages' 1 { [8] [1:3] [1:6] [1] } 1 'CCN' 0 1.2 }
%%%       'bars_avg_cfrac_strnp_TALL_CO.jpg'
%%%       }
%%%
%%%       % Stratiform Cloud fraction averages, all time points, CCN only
%%%       {
%%%       'Cloud Fraction Avg TALL CCN only, grouped by SST'
%%%       'bgraph_cfrac.h5'
%%%       { 'd' 'ST' }
%%%       'N_a (# cm^-^3)'
%%%       'Cloud Fraction'
%%%       { 'blue' 'cyan' 'magenta' }
%%%       { { 'S293', 'S298', 'S303' } 'NorthWest' }
%%%       'grouped'
%%%       { 'Averages' 1 { [12] [1:3] [1:6] [1] } 1 'CCN' 0 1.2 }
%%%       'bars_avg_cfrac_strat_TALL_CO.jpg'
%%%       }
%%%
%%%       % Convective Cloud fraction averages, all time points, CCN only
%%%       {
%%%       'Convective Cloud Fraction Avg TALL CCN only, grouped by SST'
%%%       'bgraph_cfrac.h5'
%%%       { 'a' 'CV' }
%%%       'N_a (# cm^-^3)'
%%%       'Cloud Fraction'
%%%       { 'blue' 'cyan' 'magenta' }
%%%       { { 'S293', 'S298', 'S303' } 'NorthWest' }
%%%       'grouped'
%%%       { 'Averages' 1 { [16] [1:3] [1:6] [1] } 1 'CCN' 0 1.2 }
%%%       'bars_avg_cfrac_cumul_TALL_CO.jpg'
%%%       }
%%%
%%%       % Stratiform NP-P mix Cloud fraction averages, all time points, CCN only
%%%       {
%%%       'Stratiform NP-P mix Cloud Fraction Avg TALL CCN only, grouped by SST'
%%%       'bgraph_cfrac.h5'
%%%       { 'a' 'STM' }
%%%       'N_a (# cm^-^3)'
%%%       'Cloud Fraction'
%%%       { 'blue' 'cyan' 'magenta' }
%%%       { { 'S293', 'S298', 'S303' } 'NorthWest' }
%%%       'grouped'
%%%       { 'Averages' 1 { [20] [1:3] [1:6] [1] } 1 'CCN' 0 1.2 }
%%%       'bars_avg_cfrac_stmix_TALL_CO.jpg'
%%%       }
%%%
%%%       % Stratiform-Convective mix Clouds Cloud fraction averages, all time points, CCN only
%%%       {
%%%       'All Clouds Cloud Fraction Avg TALL CCN only, grouped by SST'
%%%       'bgraph_cfrac.h5'
%%%       { 'a' 'SC' }
%%%       'N_a (# cm^-^3)'
%%%       'Cloud Fraction'
%%%       { 'blue' 'cyan' 'magenta' }
%%%       { { 'S293', 'S298', 'S303' } 'NorthWest' }
%%%       'grouped'
%%%       { 'Averages' 1 { [24] [1:3] [1:6] [1] } 1 'CCN' 0 1.2 }
%%%       'bars_avg_cfrac_scmix_TALL_CO.jpg'
%%%       }

     };


  for ipd = 1:length(PlotDefs)
    PdName  = PlotDefs{ipd}{1};
    InFile  = sprintf('%s/%s', Ddir, PlotDefs{ipd}{2});

    % Title
    Pmark   = PlotDefs{ipd}{3}{1};
    Ptitle  = PlotDefs{ipd}{3}{2};

    Xlabel  = PlotDefs{ipd}{4};
    Ylabel  = PlotDefs{ipd}{5};
    Bcolors = PlotDefs{ipd}{6};

    % LegSpec
    LegText = PlotDefs{ipd}{7}{1};
    LegLoc  = PlotDefs{ipd}{7}{2};

    Pstyle  = PlotDefs{ipd}{8};

    % VarSpec
    InVarName  = PlotDefs{ipd}{9}{1};
    InVarScale = PlotDefs{ipd}{9}{2};
    Vsel       = PlotDefs{ipd}{9}{3}{1};
    Ssel       = PlotDefs{ipd}{9}{3}{2};
    Csel       = PlotDefs{ipd}{9}{3}{3};
    Gsel       = PlotDefs{ipd}{9}{3}{4};
    Tpose      = PlotDefs{ipd}{9}{4};
    Vname      = PlotDefs{ipd}{9}{5};
    Vmin       = PlotDefs{ipd}{9}{6};
    Vmax       = PlotDefs{ipd}{9}{7};

    OutFname = PlotDefs{ipd}{10};

    fprintf('********************************************************************\n');
    fprintf('Generating bar graph: %s\n', PdName);
    fprintf('  Input File: %s\n', InFile);
    fprintf('    Variable: %s\n', InVarName);
    fprintf('    Scale: %f4\n', InVarScale);
    fprintf('  Selection specs:\n');
    fprintf('    Vsel: %s\n', strtrim(sprintf('%d ', Vsel)));
    fprintf('    Ssel: %s\n', strtrim(sprintf('%d ', Ssel)));
    fprintf('    Csel: %s\n', strtrim(sprintf('%d ', Csel)));
    fprintf('    Gsel: %s\n', strtrim(sprintf('%d ', Gsel)));
    fprintf('  Transpose: %d\n', Tpose);
    fprintf('  Vname: %s\n', Vname);
    fprintf('    Vmin: %f\n', Vmin);
    fprintf('    Vmax: %f\n', Vmax);
    fprintf('\n');
 

    % read in variables
    HDATA = hdf5read(InFile, InVarName) .* InVarScale;
    BDATA = squeeze(HDATA(Vsel, Ssel, Csel, Gsel));
    if (Tpose == 1)
      BDATA = BDATA';
    end

    XVAR  = hdf5read(InFile, Vname);

    % If doing cloud distribution:
    %   Sum up the first three columns (SNP, STM, ST)
    %   Add a count of clear space so that all bars stack up to 100%.
    % NOTE: this assumes that you've specified this plot to create percent values
    % in BDATA.
    if (regexp(PdName, '^Cloud Distribution'))
      BDATA = [ sum(BDATA(:,1:3),2) BDATA(:,4:5) ];
      % Sum up across variables (v) and subtract from 100 to get the clear
      % column counts. Then append the clear counts to the end of BDATA.
      CLEAR = 100 - nansum(BDATA, 2);
      BDATA = [ BDATA CLEAR ];
    end

    % do the plot
    clear AxisProps;
    iaxis = 0; % set this to next availble slot in the AxisProps array

    iaxis = iaxis + 1;
    AxisProps(iaxis).Name = 'FontSize';
    AxisProps(iaxis).Val  = 20;

    % x-axis labeling
    if (strcmp(Vname, 'CCN') | strcmp(Vname, 'SST'))
      % label tick marks with CCN or SST value
      Nvar = length(XVAR);
      X = 1:Nvar;      % use these for x-axis values -> evenly spaced integers creates evenly spaced bars
      clear XTIckLabels;
      for i = 1:Nvar
        XTickLabels{i} = sprintf('%d', XVAR(i));
      end
      iaxis = iaxis + 1;
      AxisProps(iaxis).Name = 'XTickLabel';
      AxisProps(iaxis).Val = XTickLabels;
    else
      X = XVAR;
    end


    iaxis = iaxis + 1;
    AxisProps(iaxis).Name = 'Ylim';
    AxisProps(iaxis).Val  = [ Vmin Vmax ];

    % plot averages
    OutFile = sprintf('%s/%s', Pdir, OutFname);
    fprintf('      %s\n', OutFile);

    PlotBarSet(X, BDATA, Ptitle, Pmark, Xlabel, Ylabel, Bcolors, Pstyle, LegText, LegLoc, AxisProps, OutFile);

    fprintf('\n');

  end
end
