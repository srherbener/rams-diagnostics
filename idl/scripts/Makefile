#
# Config
#

DEST = $(HOME)/etc/idl
SCRIPTS = $(DEST)/do_plots.pro \
          $(DEST)/plot_intensity.pro \
          $(DEST)/plot_max.pro \
          $(DEST)/plot_h_integral.pro \
          $(DEST)/plot_utils.pro \
          $(DEST)/read_revu_hdf5_var.pro \
          $(DEST)/read_grads.pro

#
# Targets
#

all: $(SCRIPTS)

$(DEST)/do_plots.pro: do_plots.pro
	cp do_plots.pro $(DEST)/do_plots.pro

$(DEST)/plot_intensity.pro: plot_intensity.pro
	cp plot_intensity.pro $(DEST)/plot_intensity.pro

$(DEST)/plot_max.pro: plot_max.pro
	cp plot_max.pro $(DEST)/plot_max.pro

$(DEST)/plot_h_integral.pro: plot_h_integral.pro
	cp plot_h_integral.pro $(DEST)/plot_h_integral.pro

$(DEST)/plot_utils.pro: plot_utils.pro
	cp plot_utils.pro $(DEST)/plot_utils.pro

$(DEST)/read_revu_hdf5_var.pro: read_revu_hdf5_var.pro
	cp read_revu_hdf5_var.pro $(DEST)/read_revu_hdf5_var.pro

$(DEST)/read_grads.pro: read_grads.pro
	cp read_grads.pro $(DEST)/read_grads.pro
