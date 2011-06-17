
SUBDIRS=src

# Where the data files are installed, generally /usr/local/share/stamp/...
stampdir=$(datadir)/@PACKAGE@
dist_stamp_DATA=defs/stamp.defaults defs/pdb.directories defs/dssp.directories

# Also distribute Tests in ./t but don't install them
EXTRA_DIST = t
dist-hook:
	rm -rf `find $(distdir)/t -name .svn`

test:
	./t/TEST
