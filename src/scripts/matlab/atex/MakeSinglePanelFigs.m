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
%  'CloudStruct_CG_C1600_S298'
%  'CloudStruct_CG_C400_G10m0'
%  'CloudStruct_CG_C400_G10m2'
%  'CloudStruct_CG_C400_G10m4'
%  'CloudStruct_CG_C400_S298'
%  'CloudStruct_CG_C400_S303'
%  'CloudStruct_CG_C50_S298'
%  'CloudStruct_CG_S298_G10M0'
%  'CloudStruct_CG_S298_G10M2'
%  'CloudStruct_CG_S298_G10M4'
%  'CloudStruct_CO_C0050_G10M5'
%  'CloudStruct_CO_C0100_G10M5'
%  'CloudStruct_CO_C0200_G10M5'
%  'CloudStruct_CO_C0400_G10M5'
%  'CloudStruct_CO_C0800_G10M5'
%  'CloudStruct_CO_C1200_G10M5'
%  'CloudStruct_CO_C1600_G10M5'
%  'CloudStruct_CO_S293_G10M5'
%  'CloudStruct_CO_S298_G10M5'
%  'CloudStruct_CO_S303_G10M5'
%  'LatHeatAvg_LWP0p01_CG_C1600_S298'
%  'LatHeatAvg_LWP0p01_CG_C400_G10M0'
%  'LatHeatAvg_LWP0p01_CG_C400_G10M2'
%  'LatHeatAvg_LWP0p01_CG_C400_G10M4'
%  'LatHeatAvg_LWP0p01_CG_C400_S298'
%  'LatHeatAvg_LWP0p01_CG_C400_S303'
%  'LatHeatAvg_LWP0p01_CG_C50_S298'
%  'LatHeatAvg_LWP0p01_CG_S298_G10M0'
%  'LatHeatAvg_LWP0p01_CG_S298_G10M2'
%  'LatHeatAvg_LWP0p01_CG_S298_G10M4'
%  'LatHeatAvg_LWP0p01_CO_C0050_G10M5'
%  'LatHeatAvg_LWP0p01_CO_C0100_G10M5'
%  'LatHeatAvg_LWP0p01_CO_C0200_G10M5'
%  'LatHeatAvg_LWP0p01_CO_C0400_G10M5'
%  'LatHeatAvg_LWP0p01_CO_C0800_G10M5'
%  'LatHeatAvg_LWP0p01_CO_C1200_G10M5'
%  'LatHeatAvg_LWP0p01_CO_C1600_G10M5'
%  'LatHeatAvg_LWP0p01_CO_S293_G10M5'
%  'LatHeatAvg_LWP0p01_CO_S298_G10M5'
%  'LatHeatAvg_LWP0p01_CO_S303_G10M5'
%  'LatHeatAvg_LWP0p10_CG_C1600_S298'
%  'LatHeatAvg_LWP0p10_CG_C400_G10M0'
%  'LatHeatAvg_LWP0p10_CG_C400_G10M2'
%  'LatHeatAvg_LWP0p10_CG_C400_G10M4'
%  'LatHeatAvg_LWP0p10_CG_C400_S298'
%  'LatHeatAvg_LWP0p10_CG_C400_S303'
%  'LatHeatAvg_LWP0p10_CG_C50_S298'
%  'LatHeatAvg_LWP0p10_CG_S298_G10M0'
%  'LatHeatAvg_LWP0p10_CG_S298_G10M2'
%  'LatHeatAvg_LWP0p10_CG_S298_G10M4'
%  'LatHeatAvg_LWP0p10_CO_C0050_G10M5'
%  'LatHeatAvg_LWP0p10_CO_C0100_G10M5'
%  'LatHeatAvg_LWP0p10_CO_C0200_G10M5'
%  'LatHeatAvg_LWP0p10_CO_C0400_G10M5'
%  'LatHeatAvg_LWP0p10_CO_C0800_G10M5'
%  'LatHeatAvg_LWP0p10_CO_C1200_G10M5'
%  'LatHeatAvg_LWP0p10_CO_C1600_G10M5'
%  'LatHeatAvg_LWP0p10_CO_S293_G10M5'
%  'LatHeatAvg_LWP0p10_CO_S298_G10M5'
%  'LatHeatAvg_LWP0p10_CO_S303_G10M5'
%  'LatHeatAvg_LWP1p00_CG_C1600_S298'
%  'LatHeatAvg_LWP1p00_CG_C400_G10M0'
%  'LatHeatAvg_LWP1p00_CG_C400_G10M2'
%  'LatHeatAvg_LWP1p00_CG_C400_G10M4'
%  'LatHeatAvg_LWP1p00_CG_C400_S298'
%  'LatHeatAvg_LWP1p00_CG_C400_S303'
%  'LatHeatAvg_LWP1p00_CG_C50_S298'
%  'LatHeatAvg_LWP1p00_CG_S298_G10M0'
%  'LatHeatAvg_LWP1p00_CG_S298_G10M2'
%  'LatHeatAvg_LWP1p00_CG_S298_G10M4'
%  'LatHeatAvg_LWP1p00_CO_C0050_G10M5'
%  'LatHeatAvg_LWP1p00_CO_C0100_G10M5'
%  'LatHeatAvg_LWP1p00_CO_C0200_G10M5'
%  'LatHeatAvg_LWP1p00_CO_C0400_G10M5'
%  'LatHeatAvg_LWP1p00_CO_C0800_G10M5'
%  'LatHeatAvg_LWP1p00_CO_C1200_G10M5'
%  'LatHeatAvg_LWP1p00_CO_C1600_G10M5'
%  'LatHeatAvg_LWP1p00_CO_S293_G10M5'
%  'LatHeatAvg_LWP1p00_CO_S298_G10M5'
%  'LatHeatAvg_LWP1p00_CO_S303_G10M5'
%  'RrateDist_CG_C1600_S298'
%  'RrateDist_CG_C400_G10M0'
%  'RrateDist_CG_C400_G10M2'
%  'RrateDist_CG_C400_G10M4'
%  'RrateDist_CG_C400_S298'
%  'RrateDist_CG_C400_S303'
%  'RrateDist_CG_C50_S298'
%  'RrateDist_CG_S298_G10M0'
%  'RrateDist_CG_S298_G10M2'
%  'RrateDist_CG_S298_G10M4'
%  'RrateDist_CO_C0050_G10M5'
%  'RrateDist_CO_C0100_G10M5'
%  'RrateDist_CO_C0200_G10M5'
%  'RrateDist_CO_C0400_G10M5'
%  'RrateDist_CO_C0800_G10M5'
%  'RrateDist_CO_C1200_G10M5'
%  'RrateDist_CO_C1600_G10M5'
%  'RrateDist_CO_S293_G10M5'
%  'RrateDist_CO_S298_G10M5'
%  'RrateDist_CO_S303_G10M5'
%  'TsAvgR2C_CG_C0400_G10M0'
%  'TsAvgR2C_CG_C1600_S298'
%  'TsAvgR2C_CG_C400_G10M2'
%  'TsAvgR2C_CG_C400_G10M4'
%  'TsAvgR2C_CG_C400_S298'
%  'TsAvgR2C_CG_C400_S303'
%  'TsAvgR2C_CG_C50_S298'
%  'TsAvgR2C_CG_S298_G10M0'
%  'TsAvgR2C_CG_S298_G10M2'
%  'TsAvgR2C_CG_S298_G10M4'
%  'TsAvgR2C_CO_C0050_G10M5'
%  'TsAvgR2C_CO_C0100_G10M5'
%  'TsAvgR2C_CO_C0200_G10M5'
%  'TsAvgR2C_CO_C0400_G10M5'
%  'TsAvgR2C_CO_C0800_G10M5'
%  'TsAvgR2C_CO_C1200_G10M5'
%  'TsAvgR2C_CO_C1600_G10M5'
%  'TsAvgR2C_CO_S293_G10M5'
%  'TsAvgR2C_CO_S298_G10M5'
%  'TsAvgR2C_CO_S303_G10M5'
%  'TsAvgTopSwup_CG_C1600_S298'
%  'TsAvgTopSwup_CG_C400_G10M0'
%  'TsAvgTopSwup_CG_C400_G10M2'
%  'TsAvgTopSwup_CG_C400_G10M4'
%  'TsAvgTopSwup_CG_C400_S298'
%  'TsAvgTopSwup_CG_C400_S303'
%  'TsAvgTopSwup_CG_C50_S298'
%  'TsAvgTopSwup_CG_S298_G10M0'
%  'TsAvgTopSwup_CG_S298_G10M2'
%  'TsAvgTopSwup_CG_S298_G10M4'
%  'TsAvgTopSwup_CO_C0050_G10M5'
%  'TsAvgTopSwup_CO_C0100_G10M5'
%  'TsAvgTopSwup_CO_C0200_G10M5'
%  'TsAvgTopSwup_CO_C0400_G10M5'
%  'TsAvgTopSwup_CO_C0800_G10M5'
%  'TsAvgTopSwup_CO_C1200_G10M5'
%  'TsAvgTopSwup_CO_C1600_G10M5'
%  'TsAvgTopSwup_CO_S293_G10M5'
%  'TsAvgTopSwup_CO_S298_G10M5'
%  'TsAvgTopSwup_CO_S303_G10M5'
%  'pop_CO_G10M5'
%  };

FigList = {
  'TsAvg_LTSS_SAMPLE'
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
