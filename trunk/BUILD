#!/bin/csh
#
#  csh script to compile the elements of STAMP
#  
#  type ./BUILD <machine>
#
#  e.g.  ./BUILD sun  or ./BUILD sgi 
#
#  This script will create a subdirectory /bin if it does not 
#  already exist.  It will then create a subsub directory /bin/machine
#  where executables for the architecture will be stored.
#
#  You should then add this directory to your path, or add links
#  from a directory that is on your path (e.g. /usr/local/bin) 
# 
echo ""
echo "Trying to BUILD STAMP for $1"
echo ""
if(! -e bin) then
	echo "Creating bin directory"
	mkdir bin
endif
echo makefile.$1
#if(! -e makefile.$1) then
#	echo "Error: no suitable makefile for "$1
#	echo "  try editing one of the src/makefile.* files yourself"
#	echo 
#	exit
#endif
if(! -e bin/$1) then
	echo  "Creating subdirectory bin/$1"
	mkdir bin/$1
endif
cd src
if(-e stamp.o) then
	echo "Deleting old object files"
	/bin/rm *.o
endif
echo "Copying Makefile"
/bin/cp makefile.$1 makefile
echo "Attempting make"
make
echo ""
echo "Moving executables to bin directory"
echo ""
/bin/mv avestruc ../bin/$1
/bin/mv alignfit ../bin/$1
/bin/mv dstamp ../bin/$1
/bin/mv extrans ../bin/$1
/bin/mv gstamp ../bin/$1
/bin/mv mergetrans ../bin/$1
/bin/mv mergestamp ../bin/$1
/bin/mv pdbc ../bin/$1
/bin/mv pdbseq ../bin/$1
/bin/mv pickframe ../bin/$1
/bin/mv poststamp ../bin/$1
/bin/mv sorttrans ../bin/$1
/bin/mv stamp ../bin/$1
/bin/mv stamp_clean ../bin/$1
/bin/mv transform ../bin/$1
/bin/mv ver2hor ../bin/$1

echo "Deleting object files"
/bin/rm  *.o

echo "Trying to put PERL programs in the right place "
which perl | awk '{ printf("\#\!%s\n",$1); }' >  ../bin/$1/aconvert
tail +2 ../perl/aconvert >> ../bin/$1/aconvert
chmod ugo+x ../bin/$1/aconvert

echo "All should be complete."
