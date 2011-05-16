#include <stdio.h>
#include <stdlib.h>
#include <stamp.h>


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
	if((DIAG & patha[j][i]) == DIAG) {
	    ++*k;
	    alseq[*k] = seqa[i-1];
	    blseq[*k] = seqb[j-1];
	    i--; j--;
	} else if((HORIZ & patha[j][i]) == HORIZ) {
	    ++*k;
	    alseq[*k] = GAP;
	    *hgap += 1;
	    blseq[*k] = seqb[j-1];
	    j--;
	} else if((VERT & patha[j][i]) == VERT) {
	    ++*k;
	    alseq[*k] = seqa[i-1];
	    blseq[*k] = GAP;
	    *vgap += 1;
	    i--;
	} else {
	    printf("Disaster in aliseq\n");
	    return 0;
	}
    }
    reval(alseq,0,*k);
    reval(blseq,0,*k);
    ++(*k);
    alseq[*k] = '\0';
    blseq[*k] = '\0';
    return 0;
}

