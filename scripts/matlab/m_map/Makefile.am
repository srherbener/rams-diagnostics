matmmapdir = $(prefix)/etc/matlab
matmmap_SCRIPTS = m_map

m_map: m_map1.4.1.tar
	tar -xzvf $<

# override install rule since we need a recursive copy, ie use
# cp -r instead of install since install does not do recursive copies
install-exec-local:
	cp -r $(matmmap_SCRIPTS) $(matmmapdir)
