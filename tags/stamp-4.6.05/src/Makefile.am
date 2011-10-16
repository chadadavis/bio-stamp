# Add this option so that data dir can be determined
AM_CFLAGS=-DAM_STAMPDIR='"$(datadir)/stamp"' -DAM_VERSION='"4.6.05"'

AM_LDFLAGS=-lm

# Build targets
bin_PROGRAMS=alignfit avestruc cofm dstamp extrans gstamp mergestamp mergetrans pdbc pdbseq pickframe poststamp sorttrans stamp stamp_clean transform ver2hor

# headers dependencies, for all binaries
noinst_HEADERS=aadat.h alignfit.h array.h cluster.h defaults.h domain.h dstamp.h f2c.h general.h gjnoc.h gjutil.h igetcb.h mergestamp.h poststamp.h rbr_aadat.h rbr_newlib.h rbr_vector.h R.h routines.h stamp.h stamprel.h ver2hor.h

# Specific deps
alignfit_SOURCES=alignfit.c a3to1.c agetbloc.c alignfitgetpars.c closefile.c count_domain.c fmatfit.c fmatmult.c get_clust.c getca.c getdomain.c getfile.c gjnoc.c gjutil.c ltou.c newoutput.c openfile.c qkfit.c readtree.c rmsp.c testfile.c update.c utol.c
avestruc_SOURCES=avestruc.c a3to1.c agetbloc.c closefile.c count_domain.c disp.c getdomain.c getfile.c getstampdat.c getstamprel.c gjutil.c igetca.c igetcb.c igetgen.c igetside.c ltou.c matmult.c openfile.c printdomain.c printmat.c rossmann.c testfile.c utol.c vector.c
cofm_SOURCES=cofm.c R.c a3to1.c closefile.c count_domain.c getdomain.c getfile.c igetca.c ltou.c matmult.c openfile.c printmat.c testfile.c utol.c
dstamp_SOURCES=dstamp.c agetbloc.c getstampdat.c getstamprel.c gjutil.c ltou.c threestate.c utol.c
extrans_SOURCES=extrans.c closefile.c count_domain.c getdomain.c getfile.c ltou.c matinv.c matprod.c matvecprod.c openfile.c printdomain.c printmat.c testfile.c utol.c
gstamp_SOURCES=gstamp.c a3to1.c agetbloc.c closefile.c count_domain.c disp.c getdomain.c getfile.c getstampdat.c getstamprel.c gjutil.c igetca.c igetcb.c igetgen.c ltou.c matmult.c openfile.c printdomain.c printmat.c rossmann.c smoothsec.c testfile.c threestate.c utol.c vector.c
mergestamp_SOURCES=mergestamp.c agetbloc.c count_domain.c getdomain.c getfile.c getstampdat.c getstamprel.c gjutil.c ltou.c matinv.c matprod.c matvecprod.c printdomain.c printmat.c rmsp.c testfile.c utol.c
mergetrans_SOURCES=mergetrans.c closefile.c count_domain.c getdomain.c getfile.c ltou.c matinv.c matprod.c matvecprod.c openfile.c printdomain.c printmat.c testfile.c utol.c
pdbc_SOURCES=pdbc.c closefile.c getfile.c ltou.c openfile.c testfile.c utol.c
pdbseq_SOURCES=pdbseq.c a3to1.c closefile.c count_domain.c getdomain.c getfile.c igetca.c ltou.c openfile.c testfile.c utol.c
pickframe_SOURCES=pickframe.c closefile.c count_domain.c getdomain.c getfile.c ltou.c matinv.c matprod.c matvecprod.c openfile.c printdomain.c printmat.c testfile.c utol.c
poststamp_SOURCES=poststamp.c a3to1.c agetbloc.c closefile.c count_domain.c disp.c getdomain.c getfile.c getstampdat.c gjutil.c igetca.c ltou.c matmult.c openfile.c printdomain.c printmat.c rossmann.c testfile.c utol.c
sorttrans_SOURCES=sorttrans.c count_domain.c getdomain.c getfile.c ltou.c printdomain.c testfile.c utol.c
stamp_clean_SOURCES=stamp_clean.c agetbloc.c clean_block.c count_domain.c getdomain.c getfile.c getstampdat.c gjutil.c ltou.c printdomain.c rmsp.c stamp_clean_block.c testfile.c utol.c
stamp_SOURCES=stamp.c a3to1.c aliseq.c ccprobcalc.c closefile.c count_domain.c display_align.c dosort.c get_clust.c getdomain.c getdsspsum.c getfile.c getks.c getsec.c gjnoc.c gjutil.c idist.c igetca.c igetcadssp.c ltou.c makeafile.c matfit.c matinv.c matmult.c murzin_P.c newoutput.c openfile.c pairfit.c pairpath.c pairwise.c printmat.c probcalc.c probplot.c qkfit.c readtree.c reval.c revmatmult.c rossmann.c roughfit.c scan.c sec_content.c slow_scan.c stampgetpars.c sw7ccs.c swstruc.c testfile.c treefit.c treepath.c treewise.c update.c utol.c
transform_SOURCES=transform.c closefile.c count_domain.c extract_dssp.c extract_pdb.c fmatmult.c getdomain.c getfile.c ltou.c openfile.c testfile.c utol.c
ver2hor_SOURCES=ver2hor.c agetbloc.c getstampdat.c getstamprel.c gjutil.c ltou.c utol.c


