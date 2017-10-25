matpivmatdir = $(prefix)/etc/matlab
matpivmat_SCRIPTS = pivmat

PM_TARFILE = pivmat.tar.gz

# use pivmat/_readmefirst_.txt as a target so that the tar file will get unpacked
pivmat: pivmat/_readmefirst_.txt $(PM_TARFILE)

pivmat/_readmefirst_.txt:
	mkdir -p pivmat
	cd pivmat; tar -xzvf ../$(PM_TARFILE)

# override install rule since we need a recursive copy, ie use
# cp -r instead of install since install does not do recursive copies
install-exec-local:
	cp -r $(matpivmat_SCRIPTS) $(matpivmatdir)
