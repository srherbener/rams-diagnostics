matnctbdir = $(prefix)/etc/matlab
matnctb_SCRIPTS = nctoolbox

NCTB_TARFILE = $(srcdir)/nctoolbox-20130305.tar.gz

# use nctoolbox/README as a target to make sure the tar file gets unpacked
nctoolbox: nctoolbox/README $(NCTB_TARFILE)

nctoolbox/README:
	mkdir -p nctoolbox
	cd nctoolbox; tar -xzvf ../$(NCTB_TARFILE)

# override install rule since we need a recursive copy, ie use
# cp -r instead of install since install does not do recursive copies
install-exec-local:
	cp -r $(matnctb_SCRIPTS) $(matnctbdir)
