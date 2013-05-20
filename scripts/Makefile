#
# Config
#

G_DEST = $(HOME)/etc/grads
G_SCRIPTS = $(G_DEST)/cbarn.gs \
            $(G_DEST)/cbarline.gs \
            $(G_DEST)/dumpStats.gs \
            $(G_DEST)/plotCfad.gs \
            $(G_DEST)/plotHslice.gs \
            $(G_DEST)/plotInitFields.gs \
            $(G_DEST)/plotRadialWind.gs \
            $(G_DEST)/plotTestCcnFields.gs \
            $(G_DEST)/dumpVar2d.gs \
            $(G_DEST)/plotVslices.gs \
            $(G_DEST)/plotTseries.gs \
            $(G_DEST)/DumpSounding.gs

B_DEST = $(HOME)/bin
B_SCRIPTS = $(B_DEST)/run_azavg_grads \
            $(B_DEST)/run_diags_grads \
            $(B_DEST)/gen_grads_plots \
            $(B_DEST)/gen_grads_stats \
            $(B_DEST)/plot_grads_vslice \
            $(B_DEST)/RunGradsBatch \
            $(B_DEST)/extractSoundings

M_DEST = $(HOME)/etc/perl_grads
#M_PACKAGES = $(M_DEST)/DiagUtils.pm

#
# Targets
#

all: $(G_SCRIPTS) $(B_SCRIPTS) $(M_PACKAGES)

$(G_DEST)/cbarn.gs: cbarn.gs
	cp cbarn.gs $(G_DEST)/cbarn.gs

$(G_DEST)/cbarline.gs: cbarline.gs
	cp cbarline.gs $(G_DEST)/cbarline.gs

$(G_DEST)/dumpStats.gs: dumpStats.gs
	cp dumpStats.gs $(G_DEST)/dumpStats.gs

$(G_DEST)/plotCfad.gs: plotCfad.gs
	cp plotCfad.gs $(G_DEST)/plotCfad.gs

$(G_DEST)/plotHslice.gs: plotHslice.gs
	cp plotHslice.gs $(G_DEST)/plotHslice.gs

$(G_DEST)/plotInitFields.gs: plotInitFields.gs
	cp plotInitFields.gs $(G_DEST)/plotInitFields.gs

$(G_DEST)/plotRadialWind.gs: plotRadialWind.gs
	cp plotRadialWind.gs $(G_DEST)/plotRadialWind.gs

$(G_DEST)/plotTestCcnFields.gs: plotTestCcnFields.gs
	cp plotTestCcnFields.gs $(G_DEST)/plotTestCcnFields.gs

$(G_DEST)/dumpVar2d.gs: dumpVar2d.gs
	cp dumpVar2d.gs $(G_DEST)/dumpVar2d.gs

$(G_DEST)/plotVslices.gs: plotVslices.gs
	cp plotVslices.gs $(G_DEST)/plotVslices.gs

$(G_DEST)/plotTseries.gs: plotTseries.gs
	cp plotTseries.gs $(G_DEST)/plotTseries.gs

$(G_DEST)/DumpSounding.gs: DumpSounding.gs
	cp DumpSounding.gs $(G_DEST)/DumpSounding.gs



$(B_DEST)/run_azavg_grads: run_azavg_grads
	cp run_azavg_grads $(B_DEST)/run_azavg_grads

$(B_DEST)/run_diags_grads: run_diags_grads
	cp run_diags_grads $(B_DEST)/run_diags_grads

$(B_DEST)/gen_grads_plots: gen_grads_plots
	cp gen_grads_plots $(B_DEST)/gen_grads_plots

$(B_DEST)/gen_grads_stats: gen_grads_stats
	cp gen_grads_stats $(B_DEST)/gen_grads_stats

$(B_DEST)/plot_grads_vslice: plot_grads_vslice
	cp plot_grads_vslice $(B_DEST)/plot_grads_vslice

$(B_DEST)/extractSoundings: extractSoundings
	cp extractSoundings $(B_DEST)/extractSoundings

$(B_DEST)/RunGradsBatch: RunGradsBatch
	cp RunGradsBatch $(B_DEST)/RunGradsBatch


$(M_DEST)/DiagUtils.pm: DiagUtils.pm
	cp DiagUtils.pm $(M_DEST)/DiagUtils.pm
