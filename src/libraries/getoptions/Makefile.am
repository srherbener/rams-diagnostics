getoptionsdir = $(prefix)/etc/fortran
getoptions_LIBRARIES = libgetoptions.a

libgetoptions_a_SOURCES = getoptions.f90

GO_MOD_FILE = getoptions.mod

# Handle the module file
install-exec-local:
	mkdir -p $(getoptionsdir)
	cp $(GO_MOD_FILE) $(getoptionsdir)

uninstall-local:
	rm -f $(getoptionsdir)/$(GO_MOD_FILE)
