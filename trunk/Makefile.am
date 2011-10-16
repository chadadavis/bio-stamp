
SUBDIRS=src

# Where the data files are installed, generally /usr/local/share/stamp/...
stampdir=$(datadir)/@PACKAGE@
dist_stamp_DATA=\
defs/stamp.defaults \
defs/pdb.directories \
defs/dssp.directories

dist_bin_SCRIPTS=\
script/stamp-update-biounit \
script/stamp-update-biounit-models \
script/stamp-update-blastdb \
script/stamp-update-dssp \
script/stamp-update-pdb \
script/stamp-update-pdbdom \
script/stamp-update-pdbseq

# Also distribute Tests in ./t but don't install them
EXTRA_DIST = t
dist-hook:
	rm -rf `find $(distdir)/t -name .svn`

test:
	./t/TEST
