#SUBDIRS=src defs
SUBDIRS=src

stampdir=$(datadir)/@PACKAGE@
dist_stamp_DATA=defs/stamp.defaults defs/pdb.directories defs/dssp.directories

test:
	./t/TEST

# Don't include .svn directories in tarball
#dist-hook:
#	rm -rf `find '$(distdir)' -name .svn`

