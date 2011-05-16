/*
Copyright (1997,1998,1999,2010) Robert B. Russell & Geoffrey J. Barton

This file is part of STAMP.

STAMP is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details. A copy of the license
can be found in the LICENSE file in the STAMP installation directory.

STAMP was developed by Robert B. Russell and Geoffrey J. Barton of
current addresses:

 Prof. Robert B. Russell (RBR)                      Prof. Geoffrey J. Barton (GJB)
 Cell Networks, University of Heidelberg            College of Life Sciences
 Room 564, Bioquant                                 University of Dundee
 Im Neuenheimer Feld 267                            Dow Street
 69120 Heidelberg                                   Dundee DD1 5EH
 Germany                                            UK
                                                
 Tel: +49 6221 54 513 62                            Tel: +44 1382 385860
 Fax: +49 6221 54 514 86                            FAX: +44 1382 385764
 Email: robert.russell@bioquant.uni-heidelberg.de   E-mail g.j.barton@dundee.ac.uk
 WWW: http://www.russell.embl-heidelberg.de         WWW: http://www.compbio.dundee.ac.uk

 All use of STAMP must cite: 

 R.B. Russell and G.J. Barton, "Multiple Protein Sequence Alignment From Tertiary
  Structure Comparison: Assignment of Global and Residue Confidence Levels",
  PROTEINS: Structure, Function, and Genetics, 14:309--323 (1992).
*/
#include <stdio.h>
#include <stdlib.h>
#include "stamp.h"


/*****************************************************************************
aliseq:  given the patha array, seqa seqb and path details, return an
alignment of the sub-regions of seqa and seqb indicated by the path.
-----------------------------------------------------------------------------*/
int aliseq(char *seqa, char *seqb, struct path *apath, unsigned char **patha, 
	   char *alseq, char *blseq, int *k, int *hgap, int *vgap) {

    unsigned char DIAG = 01; /* mask for diagonal move */
    unsigned char HORIZ= 02; /*          horizontal    */
    unsigned char VERT = 04; /*          vertical      */
    unsigned char GAP  = ' ';
    
    int i = apath->end.i;
    int j = apath->end.j;
/*
    for(ii=0; ii<apath->end.i; ++ii) {
      for(jj=0; jj<apath->end.j; ++jj) 
          printf("%1d",(int)patha[jj][ii]);
      printf("\n");
    }
*/

    *k = -1; *hgap = 0; *vgap = 0;

    while ( i > apath->start.i-1 || j > apath->start.j-1 ) {
	if((DIAG & patha[j][i]) == DIAG){
	    ++*k;
	    alseq[*k] = seqa[i-1];
	    blseq[*k] = seqb[j-1];
	    i--; j--;
	}
	else if((HORIZ & patha[j][i]) == HORIZ){
	    ++*k;
	    alseq[*k] = GAP;
	    *hgap += 1;
	    blseq[*k] = seqb[j-1];
	    j--;
	}
	else if((VERT & patha[j][i]) == VERT){
	    ++*k;
	    alseq[*k] = seqa[i-1];
	    blseq[*k] = GAP;
	    *vgap += 1;
	    i--;
	}
	else {
	    printf("Disaster in aliseq\n");
	    return 0;
	    } /* End of else... */
    }
#ifdef DBGSTEVE
/* SMJS Put in to allow printing before reverse */
    alseq[*k+1] = '\0';
    blseq[*k+1] = '\0';
    printf("Before reverse\nalseq=|%s|\nblseq=|%s|\n",alseq,blseq);
#endif
/* SMJS End added */
    reval(alseq,0,*k);
    reval(blseq,0,*k);
/* SMJS Split into two statements to fix the missing terminal residues bug */
    alseq[++*k] = '\0';
    blseq[*k] = '\0';
#ifdef DBGSTEVE
    printf("After reverse\nalseq=|%s|\nblseq=|%s|\n",alseq,blseq);
#endif
    /* TPW */
    return 1;
}

