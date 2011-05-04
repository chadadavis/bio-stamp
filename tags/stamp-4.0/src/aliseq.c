/******************************************************************************
 The computer software and associated documentation called STAMP hereinafter
 referred to as the WORK which is more particularly identified and described in 
 Appendix A of the file LICENSE.  Conditions and restrictions for use of
 this package are also in this file.

 The WORK is only available to licensed institutions.

 The WORK was developed by: 
	Robert B. Russell and Geoffrey J. Barton

 Of current addresses:

 Robert B. Russell (RBR)             Geoffrey J. Barton (GJB)
 Biomolecular Modelling Laboratory   Laboratory of Molecular Biophysics
 Imperial Cancer Research Fund       The Rex Richards Building
 Lincoln's Inn Fields, P.O. Box 123  South Parks Road
 London, WC2A 3PX, U.K.              Oxford, OX1 3PG, U.K.
 Tel: +44 171 269 3583               Tel: +44 865 275368
 FAX: +44 171 269 3417               FAX: 44 865 510454
 e-mail: russell@icrf.icnet.uk       e-mail gjb@bioch.ox.ac.uk
 WWW: http://bonsai.lif.icnet.uk/    WWW: http://geoff.biop.ox.ac.uk/

 The WORK is Copyright (1995) University of Oxford
	Administrative Offices
	Wellington Square
	Oxford OX1 2JD U.K.

 All use of the WORK must cite: 
 R.B. Russell and G.J. Barton, "Multiple Protein Sequence Alignment From Tertiary
  Structure Comparison: Assignment of Global and Residue Confidence Levels",
  PROTEINS: Structure, Function, and Genetics, 14:309--323 (1992).
*****************************************************************************/
#include <stdio.h>
#include "include.h"


/*****************************************************************************
aliseq:  given the patha array, seqa seqb and path details, return an
alignment of the sub-regions of seqa and seqb indicated by the path.
-----------------------------------------------------------------------------*/
int aliseq(seqa,seqb,apath,patha,alseq,blseq,k,hgap,vgap)

char *seqa, *seqb;
struct path *apath;
char **patha;
char *alseq, *blseq;
int *k;
int *hgap, *vgap;
{
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
    reval(alseq,0,*k);
    reval(blseq,0,*k);
    alseq[++*k] = blseq[*k] = '\0';
}

