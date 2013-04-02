#
# Config
#

DEST = $(HOME)/etc/matlab/ts_debby
DIAG_DIR = /Volumes/tasman1/sherbener/data/RAMS/TsDebby

SCRIPTS = \
  $(DEST)/PlotGrids.m \
  $(DEST)/PlotMaxWind.m \
  $(DEST)/PlotPressWind.m \
  $(DEST)/PlotSLP.m \
  $(DEST)/PlotSpeedtThetae.m \
  $(DEST)/PlotTracks.m \
  $(DEST)/SetPath.m

CONFIGS = \
  $(DIAG_DIR)/SST_CONST/DiagConfig \
  $(DIAG_DIR)/SST_OBS/DiagConfig

#
# Targets
#

all: $(DEST) $(SCRIPTS) $(CONFIGS)

$(DEST):
	mkdir -p $(DEST)

$(DEST)/%.m: %.m
	cp $(<) $(@)

$(DIAG_DIR)/SST_CONST/DiagConfig: DiagConfig.SST_CONST
	cp $(<) $(@)

$(DIAG_DIR)/SST_OBS/DiagConfig: DiagConfig.SST_OBS
	cp $(<) $(@)
