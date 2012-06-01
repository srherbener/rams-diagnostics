#
# Config
#

DEST = $(HOME)/etc/matlab
SCRIPTS = $(DEST)/GenHist2d.m \
          $(DEST)/GenRrateHist.m \
          $(DEST)/PlotHist2d.m \
          $(DEST)/PlotRrateHist.m \


#
# Targets
#

all: $(SCRIPTS)

$(DEST)/GenHist2d.m: GenHist2d.m
	cp GenHist2d.m $(DEST)/GenHist2d.m

$(DEST)/GenRrateHist.m: GenRrateHist.m
	cp GenRrateHist.m $(DEST)/GenRrateHist.m

$(DEST)/PlotHist2d.m: PlotHist2d.m
	cp PlotHist2d.m $(DEST)/PlotHist2d.m

$(DEST)/PlotRrateHist.m: PlotRrateHist.m
	cp PlotRrateHist.m $(DEST)/PlotRrateHist.m

