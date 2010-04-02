/******************************************************************************
 The computer software and associated documentation called STAMP hereinafter
 referred to as the WORK which is more particularly identified and described in 
 the LICENSE.  Conditions and restrictions for use of
 this package are also in the LICENSE.

 The WORK is only available to licensed institutions.

 The WORK was developed by: 
	Robert B. Russell and Geoffrey J. Barton

 Of current addresses:

 Robert B. Russell (RBR)	            Prof. Geoffrey J. Barton (GJB)
 EMBL Heidelberg                            School of Life Sciences
 Meyerhofstrasse 1                          University of Dundee
 D-69117 Heidelberg                         Dow Street
 Germany                                    Dundee, DD1 5EH
                                          
 Tel: +49 6221 387 473                      Tel: +44 1382 345860
 FAX: +44 6221 387 517                      FAX: +44 1382 345764
 E-mail: russell@embl-heidelberg.de         E-mail geoff@compbio.dundee.ac.uk
 WWW: http://www.russell.emb-heidelberg.de  WWW: http://www.compbio.dundee.ac.uk

   The WORK is Copyright (1997,1998,1999) Robert B. Russell & Geoffrey J. Barton
	
	
	

 All use of the WORK must cite: 
 R.B. Russell and G.J. Barton, "Multiple Protein Sequence Alignment From Tertiary
  Structure Comparison: Assignment of Global and Residue Confidence Levels",
  PROTEINS: Structure, Function, and Genetics, 14:309--323 (1992).
*****************************************************************************/
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
    
    int ii,jj;
    int i = apath->end.i;
    int j = apath->end.j;
    char temp;
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
}

